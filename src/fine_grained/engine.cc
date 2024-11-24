#include "engine.h"


int compare_long(const void* a, const void* b)
{
    int64_t arg1 = *(const int64_t*)a;
    int64_t arg2 = *(const int64_t*)b;
 
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}


DiskDataReader_R::DiskDataReader_R(EmbeddingExplorationEngine *engine)
{
    engine_ = engine;
    socket_id_ = engine->socket_id_;
    num_sockets_ = engine->num_sockets_;

    uint32_t num_nodes = engine_->Rrtree_->m_stats.getNumberOfNodes();
    for (hash_set_size_ = 1; hash_set_size_ < num_nodes; hash_set_size_ <<= 1)
    {
    }
    hash_set_size_ >>= num_confused_bits;
    node_id_mask_ = hash_set_size_ - 1;

    buff_manager_ = new NumaAwareBufferManager(MAX_NUM_SOCKETS);

    int max_depth = engine_->max_depth_;
    for (int d_i = 0; d_i <= max_depth; ++d_i)
    {
        curr_generation_[d_i] = 3;

        buff_manager_->register_numa_aware_buffer(
            sizeof(uint32_t) * hash_set_size_,
            (uint8_t **)&previous_request_generation_[d_i],
            socket_id_);

        buff_manager_->register_numa_aware_buffer(
            sizeof(id_type) * hash_set_size_,
            (uint8_t **)&previous_request_node_id_[d_i],
            socket_id_);
        buff_manager_->register_numa_aware_buffer(
            sizeof(uint8_t *) * hash_set_size_,
            (uint8_t **)&previous_request_disk_buffer_[d_i],
            socket_id_);

        // * initalize for alternative queues
        buff_manager_->register_numa_aware_buffer(
            sizeof(uint32_t) * hash_set_size_,
            (uint8_t **)&previous_request_generation_alternative_[d_i],
            socket_id_);

        buff_manager_->register_numa_aware_buffer(
            sizeof(id_type) * hash_set_size_,
            (uint8_t **)&previous_request_node_id_alternative_[d_i],
            socket_id_);
        buff_manager_->register_numa_aware_buffer(
            sizeof(uint8_t *) * hash_set_size_,
            (uint8_t **)&previous_request_disk_buffer_alternative_[d_i],
            socket_id_);

        buff_manager_->register_numa_aware_buffer(
            sizeof(DiskRef) * (num_nodes_per_request + 1),
            (uint8_t **)&requested_node_data_address[d_i],
            socket_id_);

    }
    buff_manager_->allocate_numa_aware_buffers();

    for (int d_i = 0; d_i <= max_depth; ++d_i)
    {
        memset(previous_request_generation_[d_i], 0, sizeof(uint32_t) * hash_set_size_);
        memset(previous_request_generation_alternative_[d_i], 0, sizeof(uint32_t) * hash_set_size_);
    }

    is_terminated_ = false;
    thread_ = new std::thread([&]()
                              {
                                  bind_to_core(cpu_topo[socket_id_/ENGINE_PER_SOCKET][socket_id_%ENGINE_PER_SOCKET*ENGINE_PER_SOCKET]);
                                  thread_main();
                              });
}

DiskDataReader_R::~DiskDataReader_R()
{
    is_terminated_ = true;
    //printf("Waiting for the Sender thread to join...\n");
    thread_->join();
    //printf("Sender thread joint.\n");
    delete thread_;

    buff_manager_->deallocate_numa_aware_buffers();
    delete buff_manager_;
    //printf("The request sender thread is terminated.\n");
}

int DiskDataReader_R::get_next_request_batch(lli &num_requested_nodes,
                                             lli &requested_disk_data_size,
                                             lli &requested_memcpy_size,
                                             lli &curr_embedding_idx
                                             )
{
    int max_depth = engine_->max_depth_;
    GlobalEmbeddingQueueState *global_queue_states = engine_->global_queue_states_;
    CompactExtendableEmbedding **global_embedding_queues = engine_->global_embedding_queues_;
    lli *global_embedding_queue_size = engine_->global_embedding_queue_size_;
    volatile lli *global_num_ready_embeddings = engine_->global_num_ready_embeddings1_; // here only consider global_num_ready_embeddings1_

    int depth = -1;
    for (int d_i = 0; d_i <= max_depth; ++d_i)
    {
        if (global_queue_states[d_i] == PartialReady && //@@: PartialReady 代表 global_queue 已经经过 shuffle 了
            global_num_ready_embeddings[d_i] < global_embedding_queue_size[d_i])
        {
            depth = d_i;
        }
    }
    // there is no pending extendable embeddings to perform the gather operation
    // launch another attempt later
    if (depth == -1)
    {
        return -1;
    }

    engine_->global_embedding_queue_mutex_.lock();
    if (!(global_queue_states[depth] == PartialReady &&
          global_num_ready_embeddings[depth] < global_embedding_queue_size[depth]))
    {
        engine_->global_embedding_queue_mutex_.unlock();
        return -1;
    }
    engine_->global_embedding_queue_mutex_.unlock();

    num_requested_nodes = 0;
    requested_disk_data_size = 0;
    requested_memcpy_size = 0;

    lli num_ready_embeddings = global_num_ready_embeddings[depth];
    lli num_embeddings = global_embedding_queue_size[depth];

    curr_embedding_idx = num_ready_embeddings;
    CompactExtendableEmbedding *embedding = global_embedding_queues[depth] + num_ready_embeddings;

    uint32_t curr_generation = curr_generation_[depth];
    uint32_t *previous_request_generation = previous_request_generation_[depth];
    id_type *previous_request_node_id = previous_request_node_id_[depth];
    uint8_t **previous_request_disk_buffer = previous_request_disk_buffer_[depth];

    uint32_t *previous_request_generation_alternative = previous_request_generation_alternative_[depth];
    id_type *previous_request_node_id_alternative = previous_request_node_id_alternative_[depth];
    uint8_t **previous_request_disk_buffer_alternative = previous_request_disk_buffer_alternative_[depth];

    // * We need to consider duplicated reads, merge requests
    // * because it may happen frequently
    // * do some modification here on e->new_disk_data_1
    // * have some kind of reuse strategy?

    uint8_t *curr_disk_data_ptx = embedding->new_disk_data_1;

    if (curr_generation % 2 == 1)
    {

        while (curr_embedding_idx < num_embeddings &&
               num_requested_nodes + 1 < num_nodes_per_request)
        {
            id_type id1 = embedding->id1;
            if (id1 == -1)
            {
                curr_embedding_idx = num_embeddings;
                break;
            }

            lli hash = id1 & node_id_mask_;
            bool is_duplicated = previous_request_generation[hash] == curr_generation && previous_request_node_id[hash] == id1;

            bool is_duplicated_last_batch = (previous_request_generation_alternative[hash] == curr_generation - 1) && previous_request_node_id_alternative[hash] == id1;

            lli required_disk_data_buffer_size = (!is_duplicated && !is_duplicated_last_batch) ? embedding->disk_data_length_1 : 0;

            lli required_memcpy_buffer_size = (!is_duplicated && is_duplicated_last_batch) ? embedding->disk_data_length_1 : 0;

            if (requested_memcpy_size + required_memcpy_buffer_size +
                    requested_disk_data_size + required_disk_data_buffer_size >
                disk_data_size_per_request)
            {
                break;
            }

            if (!is_duplicated && is_duplicated_last_batch)
            {
                memcpy(curr_disk_data_ptx,
                       previous_request_disk_buffer_alternative[hash],
                       embedding->disk_data_length_1);
            }

            embedding->new_disk_data_1 = is_duplicated ? previous_request_disk_buffer[hash] : curr_disk_data_ptx;

            if (!is_duplicated && !is_duplicated_last_batch)
            {
                // requested_nodes[num_requested_nodes] = id1;
                requested_node_data_address[depth][num_requested_nodes] = {embedding->new_disk_data_1, id1};
                num_requested_nodes += 1;
            }

            requested_disk_data_size += required_disk_data_buffer_size;
            requested_memcpy_size += required_memcpy_buffer_size;

            curr_disk_data_ptx += required_disk_data_buffer_size +
                                  required_memcpy_buffer_size;

            // * determine whether replace the content in the hash slot
            bool replace_the_hash_slot = previous_request_generation[hash] < curr_generation;
            assert(!(is_duplicated && replace_the_hash_slot)); // if the request is duplicated, it will never fill in the hash slot
            previous_request_generation[hash] = replace_the_hash_slot ? curr_generation : previous_request_generation[hash];
            previous_request_node_id[hash] = replace_the_hash_slot ? id1 : previous_request_node_id[hash];
            previous_request_disk_buffer[hash] = replace_the_hash_slot ? embedding->new_disk_data_1 : previous_request_disk_buffer[hash];

            ++curr_embedding_idx;
            ++embedding;
        }
    }
    else
    {

        while (curr_embedding_idx < num_embeddings &&
               num_requested_nodes + 1 < num_nodes_per_request)
        {
            id_type id1 = embedding->id1;
            if (id1 == -1)
            {
                curr_embedding_idx = num_embeddings;
                break;
            }

            lli hash = id1 & node_id_mask_;
            bool is_duplicated = previous_request_generation_alternative[hash] == curr_generation && previous_request_node_id_alternative[hash] == id1;

            bool is_duplicated_last_batch = (previous_request_generation[hash] == curr_generation - 1) && previous_request_node_id[hash] == id1;

            lli required_disk_data_buffer_size = (!is_duplicated && !is_duplicated_last_batch) ? embedding->disk_data_length_1 : 0;

            lli required_memcpy_buffer_size = (!is_duplicated && is_duplicated_last_batch) ? embedding->disk_data_length_1 : 0;

            if (requested_memcpy_size + required_memcpy_buffer_size +
                    requested_disk_data_size + required_disk_data_buffer_size >
                disk_data_size_per_request)
            {
                break;
            }

            if (!is_duplicated && is_duplicated_last_batch)
            {
                memcpy(curr_disk_data_ptx,
                       previous_request_disk_buffer[hash],
                       embedding->disk_data_length_1);
            }

            embedding->new_disk_data_1 = is_duplicated ? previous_request_disk_buffer_alternative[hash] : curr_disk_data_ptx;

            if (!is_duplicated && !is_duplicated_last_batch)
            {
                // requested_nodes[num_requested_nodes] = id1;
                requested_node_data_address[depth][num_requested_nodes] = {embedding->new_disk_data_1, id1};
                num_requested_nodes += 1;
            }

            requested_disk_data_size += required_disk_data_buffer_size;
            requested_memcpy_size += required_memcpy_buffer_size;

            curr_disk_data_ptx += required_disk_data_buffer_size +
                                  required_memcpy_buffer_size;

            // * determine whether replace the content in the hash slot
            bool replace_the_hash_slot = previous_request_generation_alternative[hash] < curr_generation;
            assert(!(is_duplicated && replace_the_hash_slot)); // if the request is duplicated, it will never fill in the hash slot
            previous_request_generation_alternative[hash] = replace_the_hash_slot ? curr_generation : previous_request_generation_alternative[hash];
            previous_request_node_id_alternative[hash] = replace_the_hash_slot ? id1 : previous_request_node_id_alternative[hash];
            previous_request_disk_buffer_alternative[hash] = replace_the_hash_slot ? embedding->new_disk_data_1 : previous_request_disk_buffer_alternative[hash];

            ++curr_embedding_idx;
            ++embedding;
        }
    }

    assert(curr_embedding_idx <= num_embeddings);
    curr_generation_[depth] += (curr_embedding_idx == num_embeddings);

    return depth;
}

void DiskDataReader_R::thread_main()
{
    assert(numa_run_on_node(socket_id_/ENGINE_PER_SOCKET) == 0);

    int max_depth = engine_->max_depth_;

    GlobalEmbeddingQueueState *global_queue_states = engine_->global_queue_states_;
    CompactExtendableEmbedding **global_embedding_queues = engine_->global_embedding_queues_;
    lli *global_embedding_queue_size = engine_->global_embedding_queue_size_;
    volatile lli *global_num_ready_embeddings = engine_->global_num_ready_embeddings1_;
    int partition_id = engine_->partition_id_;


    std::cout << "DiskDataReader_R start..." << std::endl;

    while (is_terminated_ == false)
    {
        __asm volatile("pause" ::
                           : "memory");

        lli curr_embedding_idx;
        lli requested_disk_data_size;
        lli requested_memcpy_size;
        lli num_requested_nodes;

        int depth = get_next_request_batch(num_requested_nodes,
                                           requested_disk_data_size,
                                           requested_memcpy_size,
                                           curr_embedding_idx
                                           );

        if (depth != -1)
        {
            // * first perform sort
            // std::sort(requested_node_data_address[depth], 
            //         requested_node_data_address[depth] + num_requested_nodes,
            //         [&](const DiskRef &x, const DiskRef &y) { return x.node_id < y.node_id; });

            // * read data from disk
            for (lli i = 0; i < num_requested_nodes; ++i)
            {
                // id_type node_id = requested_nodes[i];
                id_type node_id = requested_node_data_address[depth][i].node_id;
                uint8_t *node_data_address = requested_node_data_address[depth][i].addr;
                // * call newly-implemented readNode()
                engine_->Rrtree_->readRawData(node_id, node_data_address);
            }

            engine_->global_embedding_queue_mutex_.lock();
            global_num_ready_embeddings[depth] = curr_embedding_idx;
            engine_->global_embedding_queue_mutex_.unlock();

            total_data_read += requested_disk_data_size;
        }
    }
    // * numa_free all arrays
    // numa_free(requested_nodes, sizeof(id_type) * (num_nodes_per_request + 1));
}

//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================

DiskDataReader_S::DiskDataReader_S(EmbeddingExplorationEngine *engine)
{
    engine_ = engine;
    socket_id_ = engine->socket_id_;
    num_sockets_ = engine->num_sockets_;

    uint32_t num_nodes = engine_->Srtree_->m_stats.getNumberOfNodes();
    for (hash_set_size_ = 1; hash_set_size_ < num_nodes; hash_set_size_ <<= 1)
    {
    }
    hash_set_size_ >>= num_confused_bits;
    node_id_mask_ = hash_set_size_ - 1;

    buff_manager_ = new NumaAwareBufferManager(MAX_NUM_SOCKETS);

    int max_depth = engine_->max_depth_;
    for (int d_i = 0; d_i <= max_depth; ++d_i)
    {
        curr_generation_[d_i] = 3;

        buff_manager_->register_numa_aware_buffer(
            sizeof(uint32_t) * hash_set_size_,
            (uint8_t **)&previous_request_generation_[d_i],
            socket_id_);
        buff_manager_->register_numa_aware_buffer(
            sizeof(id_type) * hash_set_size_,
            (uint8_t **)&previous_request_node_id_[d_i],
            socket_id_);
        buff_manager_->register_numa_aware_buffer(
            sizeof(uint8_t *) * hash_set_size_,
            (uint8_t **)&previous_request_disk_buffer_[d_i],
            socket_id_);
        // * for alternative queues
        buff_manager_->register_numa_aware_buffer(
            sizeof(uint32_t) * hash_set_size_,
            (uint8_t **)&previous_request_generation_alternative_[d_i],
            socket_id_);
        buff_manager_->register_numa_aware_buffer(
            sizeof(id_type) * hash_set_size_,
            (uint8_t **)&previous_request_node_id_alternative_[d_i],
            socket_id_);
        buff_manager_->register_numa_aware_buffer(
            sizeof(uint8_t *) * hash_set_size_,
            (uint8_t **)&previous_request_disk_buffer_alternative_[d_i],
            socket_id_);

        buff_manager_->register_numa_aware_buffer(
            sizeof(DiskRef) * (num_nodes_per_request + 1),
            (uint8_t **)&requested_node_data_address[d_i],
            socket_id_);
    }
    buff_manager_->allocate_numa_aware_buffers();

    for (int d_i = 0; d_i <= max_depth; ++d_i)
    {
        memset(previous_request_generation_[d_i], 0, sizeof(uint32_t) * hash_set_size_);
        memset(previous_request_generation_alternative_[d_i], 0, sizeof(uint32_t) * hash_set_size_);
    }

    is_terminated_ = false;
    thread_ = new std::thread([&]()
                              {
                                  bind_to_core(cpu_topo[socket_id_/ENGINE_PER_SOCKET][socket_id_%ENGINE_PER_SOCKET*ENGINE_PER_SOCKET+1]);
                                  thread_main();
                              });
}

DiskDataReader_S::~DiskDataReader_S()
{
    is_terminated_ = true;
    //printf("Waiting for the Sender thread to join...\n");
    thread_->join();
    //printf("Sender thread joint.\n");
    delete thread_;

    buff_manager_->deallocate_numa_aware_buffers();
    delete buff_manager_;
    //printf("The request sender thread is terminated.\n");
}

int DiskDataReader_S::get_next_request_batch(lli &num_requested_nodes,
                                             lli &requested_disk_data_size,
                                             lli &requested_memcpy_size,
                                             lli &curr_embedding_idx
                                             )
{
    int max_depth = engine_->max_depth_;
    GlobalEmbeddingQueueState *global_queue_states = engine_->global_queue_states_;
    CompactExtendableEmbedding **global_embedding_queues = engine_->global_embedding_queues_;
    lli *global_embedding_queue_size = engine_->global_embedding_queue_size_;
    volatile lli *global_num_ready_embeddings = engine_->global_num_ready_embeddings2_; // here only consider global_num_ready_embeddings1_

    int depth = -1;
    for (int d_i = 0; d_i <= max_depth; ++d_i)
    {
        if (global_queue_states[d_i] == PartialReady && //@@: PartialReady 代表 global_queue 已经经过 shuffle 了
            global_num_ready_embeddings[d_i] < global_embedding_queue_size[d_i])
        {
            depth = d_i;
        }
    }
    // there is no pending extendable embeddings to perform the gather operation
    // launch another attempt later
    if (depth == -1)
    {
        return -1;
    }
    engine_->global_embedding_queue_mutex_.lock();
    if (!(global_queue_states[depth] == PartialReady &&
          global_num_ready_embeddings[depth] < global_embedding_queue_size[depth]))
    {
        engine_->global_embedding_queue_mutex_.unlock();
        return -1;
    }
    engine_->global_embedding_queue_mutex_.unlock();

    num_requested_nodes = 0;
    requested_disk_data_size = 0;
    requested_memcpy_size = 0;

    lli num_ready_embeddings = global_num_ready_embeddings[depth];
    lli num_embeddings = global_embedding_queue_size[depth];
    curr_embedding_idx = num_ready_embeddings;
    CompactExtendableEmbedding *embedding = global_embedding_queues[depth] + num_ready_embeddings;

    uint32_t curr_generation = curr_generation_[depth];
    uint32_t *previous_request_generation = previous_request_generation_[depth];
    id_type *previous_request_node_id = previous_request_node_id_[depth];
    uint8_t **previous_request_disk_buffer = previous_request_disk_buffer_[depth];

    uint32_t *previous_request_generation_alternative = previous_request_generation_alternative_[depth];
    id_type *previous_request_node_id_alternative = previous_request_node_id_alternative_[depth];
    uint8_t **previous_request_disk_buffer_alternative = previous_request_disk_buffer_alternative_[depth];

    // * We need to consider duplicated reads, merge requests
    // * because it may happen frequently
    // * do some modification here on e->new_disk_data_2
    // * have some kind of reuse strategy?

    uint8_t *curr_disk_data_ptx = embedding->new_disk_data_2;

    if (curr_generation % 2 == 1)
    {

        while (curr_embedding_idx < num_embeddings &&
               num_requested_nodes + 1 < num_nodes_per_request)
        {
            id_type id2 = embedding->id2;
            if (id2 == -1)
            {
                curr_embedding_idx = num_embeddings;
                break;
            }

            lli hash = id2 & node_id_mask_;
            bool is_duplicated = previous_request_generation[hash] == curr_generation && previous_request_node_id[hash] == id2;

            bool is_duplicated_last_batch = (previous_request_generation_alternative[hash] == curr_generation - 1) && previous_request_node_id_alternative[hash] == id2;

            lli required_disk_data_buffer_size = (!is_duplicated && !is_duplicated_last_batch) ? embedding->disk_data_length_2 : 0;

            lli required_memcpy_buffer_size = (!is_duplicated && is_duplicated_last_batch) ? embedding->disk_data_length_2 : 0;

            if (requested_memcpy_size + required_memcpy_buffer_size +
                    requested_disk_data_size + required_disk_data_buffer_size >
                disk_data_size_per_request)
            {
                break;
            }

            if (!is_duplicated && is_duplicated_last_batch)
            {
                memcpy(curr_disk_data_ptx,
                       previous_request_disk_buffer_alternative[hash],
                       embedding->disk_data_length_2);
            }

            embedding->new_disk_data_2 = is_duplicated ? previous_request_disk_buffer[hash] : curr_disk_data_ptx;

            if (!is_duplicated && !is_duplicated_last_batch)
            {
                // requested_nodes[num_requested_nodes] = id2;
                requested_node_data_address[depth][num_requested_nodes] = {embedding->new_disk_data_2, id2};
                num_requested_nodes += 1;
            }

            requested_disk_data_size += required_disk_data_buffer_size;
            requested_memcpy_size += required_memcpy_buffer_size;

            curr_disk_data_ptx += required_disk_data_buffer_size +
                                  required_memcpy_buffer_size;

            // * determine whether replace the content in the hash slot
            bool replace_the_hash_slot = previous_request_generation[hash] < curr_generation;
            assert(!(is_duplicated && replace_the_hash_slot)); // if the request is duplicated, it will never fill in the hash slot
            previous_request_generation[hash] = replace_the_hash_slot ? curr_generation : previous_request_generation[hash];
            previous_request_node_id[hash] = replace_the_hash_slot ? id2 : previous_request_node_id[hash];
            previous_request_disk_buffer[hash] = replace_the_hash_slot ? embedding->new_disk_data_2 : previous_request_disk_buffer[hash];

            ++curr_embedding_idx;
            ++embedding;
        }
    }
    else
    {

        while (curr_embedding_idx < num_embeddings &&
               num_requested_nodes + 1 < num_nodes_per_request)
        {
            id_type id2 = embedding->id2;
            if (id2 == -1)
            {
                curr_embedding_idx = num_embeddings;
                break;
            }

            lli hash = id2 & node_id_mask_;
            bool is_duplicated = previous_request_generation_alternative[hash] == curr_generation && previous_request_node_id_alternative[hash] == id2;

            bool is_duplicated_last_batch = (previous_request_generation[hash] == curr_generation - 1) && previous_request_node_id[hash] == id2;

            lli required_disk_data_buffer_size = (!is_duplicated && !is_duplicated_last_batch) ? embedding->disk_data_length_2 : 0;

            lli required_memcpy_buffer_size = (!is_duplicated && is_duplicated_last_batch) ? embedding->disk_data_length_2 : 0;

            if (requested_memcpy_size + required_memcpy_buffer_size +
                    requested_disk_data_size + required_disk_data_buffer_size >
                disk_data_size_per_request)
            {
                break;
            }

            if (!is_duplicated && is_duplicated_last_batch)
            {
                memcpy(curr_disk_data_ptx,
                       previous_request_disk_buffer[hash],
                       embedding->disk_data_length_2);
            }

            embedding->new_disk_data_2 = is_duplicated ? previous_request_disk_buffer_alternative[hash] : curr_disk_data_ptx;

            if (!is_duplicated && !is_duplicated_last_batch)
            {
                // requested_nodes[num_requested_nodes] = id2;
                requested_node_data_address[depth][num_requested_nodes] = {embedding->new_disk_data_2, id2};
                num_requested_nodes += 1;
            }

            requested_disk_data_size += required_disk_data_buffer_size;
            requested_memcpy_size += required_memcpy_buffer_size;

            curr_disk_data_ptx += required_disk_data_buffer_size +
                                  required_memcpy_buffer_size;

            // * determine whether replace the content in the hash slot
            bool replace_the_hash_slot = previous_request_generation_alternative[hash] < curr_generation;
            assert(!(is_duplicated && replace_the_hash_slot)); // if the request is duplicated, it will never fill in the hash slot
            previous_request_generation_alternative[hash] = replace_the_hash_slot ? curr_generation : previous_request_generation_alternative[hash];
            previous_request_node_id_alternative[hash] = replace_the_hash_slot ? id2 : previous_request_node_id_alternative[hash];
            previous_request_disk_buffer_alternative[hash] = replace_the_hash_slot ? embedding->new_disk_data_2 : previous_request_disk_buffer_alternative[hash];

            ++curr_embedding_idx;
            ++embedding;
        }
    }

    assert(curr_embedding_idx <= num_embeddings);
    curr_generation_[depth] += (curr_embedding_idx == num_embeddings);

    return depth;
}

void DiskDataReader_S::thread_main()
{
    assert(numa_run_on_node(socket_id_/ENGINE_PER_SOCKET) == 0);

    int max_depth = engine_->max_depth_;

    GlobalEmbeddingQueueState *global_queue_states = engine_->global_queue_states_;
    CompactExtendableEmbedding **global_embedding_queues = engine_->global_embedding_queues_;
    lli *global_embedding_queue_size = engine_->global_embedding_queue_size_;
    volatile lli *global_num_ready_embeddings = engine_->global_num_ready_embeddings2_;
    int partition_id = engine_->partition_id_;

    std::cout << "DiskDataReader_S start..." << std::endl;

    while (is_terminated_ == false)
    {
        __asm volatile("pause" ::
                           : "memory");

        lli curr_embedding_idx;
        lli requested_disk_data_size;
        lli requested_memcpy_size;
        lli num_requested_nodes;

        int depth = get_next_request_batch(num_requested_nodes,
                                           requested_disk_data_size,
                                           requested_memcpy_size,
                                           curr_embedding_idx
                                           );
        if (depth != -1)
        {
            // * first perform sort
            // std::sort(requested_node_data_address[depth], 
            //         requested_node_data_address[depth] + num_requested_nodes,
            //         [&](const DiskRef &x, const DiskRef &y) { return x.node_id < y.node_id; });

            // * read data from disk
            for (lli i = 0; i < num_requested_nodes; ++i)
            {
                // id_type node_id = requested_nodes[i];
                id_type node_id = requested_node_data_address[depth][i].node_id;
                uint8_t *node_data_address = requested_node_data_address[depth][i].addr;
                // * call newly-implemented readNode()
                engine_->Srtree_->readRawData(node_id, node_data_address);
            }

            engine_->global_embedding_queue_mutex_.lock();
            global_num_ready_embeddings[depth] = curr_embedding_idx;
            engine_->global_embedding_queue_mutex_.unlock();

            total_data_read += requested_disk_data_size;
        }
    }
    // * numa_free all arrays
    // numa_free(requested_nodes, sizeof(id_type) * (num_nodes_per_request + 1));
}

//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================

EmbeddingExplorationEngine::EmbeddingExplorationEngine(
    lli max_depth,
    lli max_depth_R,
    lli max_depth_S,
    lli max_num_embeddings_global_queue,
    lli global_disk_data_buffer_size,
    int num_partitions,
    int socket_id,
    int num_sockets,
    int num_computation_threads,
    int chunk_size,
    RTreeType *Rrtree,
    RTreeType *Srtree,
    Aggregator<size_t> *count)
{
    num_sockets_ = num_sockets;
    socket_id_ = socket_id;
    chunk_size_ = chunk_size;

    num_partitions_ = num_partitions;

    Rrtree_ = Rrtree;
    Srtree_ = Srtree;

    num_computation_threads_ = num_computation_threads;
    std::cout << "      num_computation_threads_ = " << num_computation_threads_ << std::endl;
    const int SCALE = num_computation_threads_ * 8;

    max_depth_ = max_depth;
    max_depth_R_ = max_depth_R;
    max_depth_S_ = max_depth_S;

    max_num_embeddings_global_queue_ = max_num_embeddings_global_queue;
    global_disk_data_buffer_size_ = global_disk_data_buffer_size;

    max_num_embeddings_local_queue_ = L1_DCACHE_SIZE / 2 / sizeof(CompactExtendableEmbedding);
    local_disk_data_buffer_size_ = global_disk_data_buffer_size / SCALE;

    count_ = count;

    // initialize all the buffers
    init_buffers();

    // clear all the buffers
    clear_all_buffers();

    // start disk reader threads
    disk_reader1_ = new DiskDataReader_R(this);
    disk_reader2_ = new DiskDataReader_S(this);

    for (int depth = 0; depth <= max_depth_; ++depth)
    {
        pthread_barrier_init(&comp_thread_barrier_[depth], NULL, num_computation_threads_ + 1);
    }
}

EmbeddingExplorationEngine::~EmbeddingExplorationEngine()
{
    delete disk_reader1_;
    delete disk_reader2_;
    release_buffer();
}

void EmbeddingExplorationEngine::init_buffers()
{
    buff_manager_ = new NumaAwareBufferManager(MAX_NUM_SOCKETS);
    int s_i = socket_id_;

    for (int d_i = 0; d_i <= max_depth_; ++d_i)
    {
        buff_manager_->register_numa_aware_buffer(
            sizeof(CompactExtendableEmbedding) * max_num_embeddings_global_queue_,
            (uint8_t **)&global_embedding_queues_[d_i], s_i);
    }

    // make sure that the graph data buffer is a continous space
    for (int d_i = 0; d_i <= max_depth_; ++d_i)
    {
        buff_manager_->register_numa_aware_buffer(global_disk_data_buffer_size_,
                                                  &global_disk_data_buffer1_[d_i], s_i);
    }

    for (int d_i = 0; d_i <= max_depth_; ++d_i)
    {
        buff_manager_->register_numa_aware_buffer(global_disk_data_buffer_size_,
                                                  &global_disk_data_buffer2_[d_i], s_i);
    }
    // * for alternative queues
    for (int d_i = 0; d_i <= max_depth_; ++d_i)
    {
        buff_manager_->register_numa_aware_buffer(global_disk_data_buffer_size_,
                                                  &global_disk_data_buffer1_alternative_[d_i], s_i);
    }

    for (int d_i = 0; d_i <= max_depth_; ++d_i)
    {
        buff_manager_->register_numa_aware_buffer(global_disk_data_buffer_size_,
                                                  &global_disk_data_buffer2_alternative_[d_i], s_i);
    }

    for (int t_i = 0; t_i < num_computation_threads_; ++t_i)
    {
        for (int d_i = 0; d_i <= max_depth_; ++d_i)
        {
            buff_manager_->register_numa_aware_buffer(
                sizeof(CompactExtendableEmbedding) * max_num_embeddings_local_queue_,
                (uint8_t **)&local_embedding_queues_[t_i][d_i], s_i);
        }
    }

    buff_manager_->allocate_numa_aware_buffers();

    lli engine_buff_size = buff_manager_->get_total_allocated_size();

    printf("*** Each engine takes %.3f (GB) buffer size\n",
           engine_buff_size / 1024. / 1024. / 1024.);
}

void EmbeddingExplorationEngine::release_buffer()
{
    buff_manager_->deallocate_numa_aware_buffers();
    delete buff_manager_;
}

void EmbeddingExplorationEngine::clear_all_buffers()
{
    //throw_entry_executable_counter_ = 0;
    for (int i = 0; i <= max_depth_; ++i)
    {
        global_queue_states_[i] = InValid;
        clear_global_queue(i);
    }
    for (int i = 0; i < num_computation_threads_; ++i)
    {
        for (int j = 0; j <= max_depth_; ++j)
        {
            local_embedding_queue_size_[i][j] = 0;
            local_needed_disk_data_buffer_size1_[i][j] = 0;
            local_needed_disk_data_buffer_size2_[i][j] = 0;
        }
    }
}

void EmbeddingExplorationEngine::flush_all_extendable_embeddings()
{
    change_global_queue_state_to_partial_ready(0);
    extend_embeddings(0);
}

void EmbeddingExplorationEngine::change_global_queue_state_to_partial_ready(int depth)
{
    assert(global_queue_states_[depth] == InValid);
    assert(disk_reader1_->curr_generation_[depth] == disk_reader2_->curr_generation_[depth]);

    uint32_t controller = disk_reader1_->curr_generation_[depth];

    if (controller % 2 == 1)
        shuffle_global_embedding_queue(depth);
    else
        shuffle_global_embedding_alternative_queue(depth);

    global_embedding_queue_mutex_.lock();
    global_queue_states_[depth] = PartialReady;
    global_embedding_queue_mutex_.unlock();
}

void EmbeddingExplorationEngine::shuffle_global_embedding_queue(int depth)
{
    lli num_embeddings = global_embedding_queue_size_[depth];
    CompactExtendableEmbedding *embedding = global_embedding_queues_[depth];

    lli previous_disk_data_size;
    lli accumlator = 0;
    for (lli i = 0; i < num_embeddings; ++i, ++embedding)
    {
        previous_disk_data_size = embedding->disk_data_length_1; // in byte
        embedding->new_disk_data_1 = global_disk_data_buffer1_[depth] + accumlator;
        accumlator += previous_disk_data_size;
    }
    assert(accumlator <= global_disk_data_buffer_size_);

    // * reset iterators
    embedding = global_embedding_queues_[depth];
    accumlator = 0;
    for (lli i = 0; i < num_embeddings; ++i, ++embedding)
    {
        previous_disk_data_size = embedding->disk_data_length_2; // in byte
        embedding->new_disk_data_2 = global_disk_data_buffer2_[depth] + accumlator;
        accumlator += previous_disk_data_size;
    }
    assert(accumlator < global_disk_data_buffer_size_);

    global_num_ready_embeddings1_[depth] = 0;
    global_num_ready_embeddings2_[depth] = 0;
}

void EmbeddingExplorationEngine::shuffle_global_embedding_alternative_queue(int depth)
{
    lli num_embeddings = global_embedding_queue_size_[depth];
    CompactExtendableEmbedding *embedding = global_embedding_queues_[depth];

    lli previous_disk_data_size;
    lli accumlator = 0;
    for (lli i = 0; i < num_embeddings; ++i, ++embedding)
    {
        previous_disk_data_size = embedding->disk_data_length_1; // in byte
        embedding->new_disk_data_1 = global_disk_data_buffer1_alternative_[depth] + accumlator;
        accumlator += previous_disk_data_size;
    }
    assert(accumlator <= global_disk_data_buffer_size_);

    // * reset iterators
    embedding = global_embedding_queues_[depth];
    accumlator = 0;
    for (lli i = 0; i < num_embeddings; ++i, ++embedding)
    {
        previous_disk_data_size = embedding->disk_data_length_2; // in byte
        embedding->new_disk_data_2 = global_disk_data_buffer2_alternative_[depth] + accumlator;
        accumlator += previous_disk_data_size;
    }
    assert(accumlator < global_disk_data_buffer_size_);

    global_num_ready_embeddings1_[depth] = 0;
    global_num_ready_embeddings2_[depth] = 0;
}

void EmbeddingExplorationEngine::extend_embeddings(int depth)
{
    // std::cout << "############## depth = " << depth << std::endl;
    assert(depth >= 0 && depth <= max_depth_);
    assert(global_queue_states_[depth] == PartialReady);

    // next_depth_ = depth + 1;
    volatile lli num_distributed_embeddings = 0; // the number of embeddings that have been assigned to a thread
    volatile lli num_extended_embeddings = 0;    // the number of extended (processed) embeddings
    lli num_embeddings_to_extend = global_embedding_queue_size_[depth];

    volatile bool is_terminated = false;

    num_suspended_threads_[depth] = 0;
    global_phases_[depth] = 0;

    // std::cout << "HERE......." << std::endl;

    auto computation_thread_main = [&, depth](int thread_id)
    {
        assert(numa_run_on_node(socket_id_/ENGINE_PER_SOCKET) == 0);
        bind_to_core(cpu_topo[socket_id_/ENGINE_PER_SOCKET][which_core_to_use(socket_id_, thread_id)]);

        SimpleNodePtr sn1, sn2;

        lli thread_begin = 0;
        lli thread_curr = 0;
        lli thread_end = 0;

        local_phases_[depth][thread_id][0] = 0; //@@ 0: idle, 1: working

        while (true)
        {
            // suspend the computation thread
            {
                std::lock_guard<std::mutex> lk(num_suspended_comp_threads_mutex_[depth]);
                if (global_phases_[depth] != local_phases_[depth][thread_id][0])
                {
                    global_phases_[depth] = local_phases_[depth][thread_id][0];
                    num_suspended_threads_[depth] = 1;
                }
                else
                {
                    num_suspended_threads_[depth] += 1;
                }
                if (num_suspended_threads_[depth] == num_computation_threads_)
                {
                    // this indicates that all comp threads will go to sleep
                    // so that main thread should be woken up
                    num_suspended_comp_threads_cv_[depth].notify_one();
                }
            }
            pthread_barrier_wait(&comp_thread_barrier_[depth]);
            local_phases_[depth][thread_id][0] ^= 1;
            if (is_terminated)
            { // @@ a local "volatile" variable
                break;
            }

            // std::cout << "COMPER....... depth = " << depth << std::endl;

            while (true)
            {
                if (thread_curr >= thread_end)
                {
                    thread_begin = thread_curr = __sync_fetch_and_add(&num_distributed_embeddings, chunk_size_);
                    thread_end = (thread_begin + chunk_size_) < num_embeddings_to_extend ? (thread_begin + chunk_size_) : num_embeddings_to_extend;
                } //@@ add a parenthesis here is more readable
                if (thread_curr >= thread_end)
                {
                    // this indicates that there is no more pending workload
                    // the computation thread should go to sleep and wait for a new workload batch
                    break;
                }

                th_wait_time[thread_id] -= get_time();
                while (num_partitions_ >= 1)
                {
                    __asm volatile("pause" ::
                                       : "memory"); // self-loop
                    bool cond = std::min(global_num_ready_embeddings1_[depth],
                                         global_num_ready_embeddings2_[depth]) >= thread_end;
                    if (cond)
                        break;
                }
                th_wait_time[thread_id] += get_time();

                // std::cout << "ready # = " << std::min(global_num_ready_embeddings1_[depth],
                //                          global_num_ready_embeddings2_[depth]) << ", req = " << thread_end << std::endl;

                // start actual computation
                CompactExtendableEmbedding *compact_e = global_embedding_queues_[depth] + thread_curr;
                while (thread_curr < thread_end)
                {
                    // * do spatial join here
                    // * if tree height not equal, the embedding like {-1, another_new_id}

                    // NodePtr * n1 = (NodePtr *)compact_e->new_disk_data_1;
                    // NodePtr * n2 = (NodePtr *)compact_e->new_disk_data_2;
                    // spatial_join_func(compact_e->id1, compact_e->id2, *n1, *n2, thread_id, depth);

                    // data is correctly read from the disk
                    // we need a simplified "NodePtr" struct to avoid use "NodePtr" of Rtree (Errors!)

                    if (compact_e->id1 == -1 && compact_e->id2 == -1)
                    {
                        throw std::runtime_error("Error...\n");
                    }
                    else if (compact_e->id1 == -1)
                    {
                        sn2.fill_data(compact_e->id2, compact_e->new_disk_data_2);
                        spatial_join_func_simplenode_oneside(compact_e, &sn2, thread_id, depth);
                    }
                    else if (compact_e->id2 == -1)
                    {
                        sn1.fill_data(compact_e->id1, compact_e->new_disk_data_1);
                        spatial_join_func_simplenode_oneside(compact_e, &sn1, thread_id, depth);
                    }
                    else
                    {
                        sn1.fill_data(compact_e->id1, compact_e->new_disk_data_1);
                        sn2.fill_data(compact_e->id2, compact_e->new_disk_data_2);
                        spatial_join_func_simplenode(compact_e->id1, compact_e->id2, &sn1, &sn2, thread_id, depth);
                    }

                    ++thread_curr;
                    ++compact_e;
                }
                lli num_locally_extended_embeddings = thread_end - thread_begin;
                __sync_fetch_and_add(&num_extended_embeddings, num_locally_extended_embeddings);
            }
        }
    };

    std::thread *computation_threads[num_computation_threads_];
    for (int t_i = 0; t_i < num_computation_threads_; ++t_i)
    {
        computation_threads[t_i] = new std::thread(computation_thread_main, t_i);
    }

    // std::cout << "HERE....... depth = " << depth << std::endl;

    int main_thread_local_phase = 0;
    while (num_extended_embeddings < num_embeddings_to_extend)
    {

        // std::cout << "HERE....... depth = " << depth << std::endl;

        pthread_barrier_wait(&comp_thread_barrier_[depth]);
        main_thread_local_phase ^= 1;
        // wait until all computation threads are suspended
        std::unique_lock<std::mutex> lk(num_suspended_comp_threads_mutex_[depth]);
        num_suspended_comp_threads_cv_[depth].wait(
            lk, [&]
            { return global_phases_[depth] == main_thread_local_phase &&
                     num_suspended_threads_[depth] == num_computation_threads_; });
        lk.unlock();

        // std::cout << "HERE....... depth = " << depth << std::endl;

        if (depth < max_depth_)
        {
            // since some embeddings are not completely extended
            // it falls within the first case
            // std::cout << "HERE....... depth = " << depth << std::endl;
            if (num_extended_embeddings < num_embeddings_to_extend)
            {
                // the next-level should not be empty
                assert(global_embedding_queue_size_[depth + 1] > 0);
                change_global_queue_state_to_partial_ready(depth + 1);
                // std::cout << "HERE....... depth = " << depth << std::endl;
                extend_embeddings(depth + 1);

                // next_depth_ = depth + 1;
                // std::cout << "next_depth is now " << depth + 1 << std::endl;
                // at this point, the next-level queue should be empty
                assert(global_embedding_queue_size_[depth + 1] == 0);
            }
        }
        else
        {
            // if the current level is the last level
            // the first case should not occur
            assert(num_extended_embeddings == num_embeddings_to_extend);
        }
    }

    // so far, all embeddings at this level have been processed
    // wake up the computation threads and terminate them
    is_terminated = true;
    pthread_barrier_wait(&comp_thread_barrier_[depth]);
    for (int t_i = 0; t_i < num_computation_threads_; ++t_i)
    {
        computation_threads[t_i]->join();
        delete computation_threads[t_i];
    }

    // at this point, the next-level queue may not be empty
    // we should flush it
    if (depth < max_depth_)
    {
        while (true)
        {
            // flush the local buffer first
            for (int t_i = 0; t_i < num_computation_threads_; ++t_i)
            {
                flush_local_embeddings(t_i, depth + 1);
            }
            if (global_embedding_queue_size_[depth + 1] > 0)
            {
                change_global_queue_state_to_partial_ready(depth + 1);
                extend_embeddings(depth + 1);
                // next_depth_ = depth + 1;
                // std::cout << "next_depth is now " << next_depth_ << std::endl;
            }
            else
            {
                break;
            }
        }
    }

    // there should be no more embeddings in the next-level queue
    // clear to current-level embedding queue
    global_queue_states_[depth] = InValid;
    clear_global_queue(depth);

    // std::cout << "Extend_embedding() is done... " << std::endl;
}

int EmbeddingExplorationEngine::flush_local_embeddings(int thread_id, int depth)
{
    if (depth <= max_depth_ && local_embedding_queue_size_[thread_id][depth] > 0)
    {

        lli local_queue_size = local_embedding_queue_size_[thread_id][depth];
        lli local_used_disk_data_size1 = local_needed_disk_data_buffer_size1_[thread_id][depth];
        lli local_used_disk_data_size2 = local_needed_disk_data_buffer_size2_[thread_id][depth];

        lli new_embedding_queue_size;
        lli new_used_disk_data_size1;
        lli new_used_disk_data_size2;

        global_embedding_queue_mutex_.lock();
        new_embedding_queue_size = global_embedding_queue_size_[depth] + local_queue_size;
        new_used_disk_data_size1 = global_used_disk_data_buffer_size1_[depth] + local_used_disk_data_size1;
        new_used_disk_data_size2 = global_used_disk_data_buffer_size2_[depth] + local_used_disk_data_size2;

        if (new_embedding_queue_size <= max_num_embeddings_global_queue_ &&
            new_used_disk_data_size1 <= global_disk_data_buffer_size_ &&
            new_used_disk_data_size2 <= global_disk_data_buffer_size_)
        {
            global_embedding_queue_size_[depth] = new_embedding_queue_size;
            global_used_disk_data_buffer_size1_[depth] = new_used_disk_data_size1;
            global_used_disk_data_buffer_size2_[depth] = new_used_disk_data_size2;
            global_embedding_queue_mutex_.unlock();
        }
        else
        {
            global_embedding_queue_mutex_.unlock();
            bool exceed_queue_size = new_embedding_queue_size > max_num_embeddings_global_queue_;
            bool exceed_disk_buffer_size1 = new_used_disk_data_size1 > global_disk_data_buffer_size_;
            bool exceed_disk_buffer_size2 = new_used_disk_data_size2 > global_disk_data_buffer_size_;
            return -1;
        }

        // copying local data to the global queue
        memcpy(&global_embedding_queues_[depth][new_embedding_queue_size - local_queue_size],
               &local_embedding_queues_[thread_id][depth][0],
               sizeof(CompactExtendableEmbedding) * local_queue_size);

        local_embedding_queue_size_[thread_id][depth] = 0;
        local_needed_disk_data_buffer_size1_[thread_id][depth] = 0;
        local_needed_disk_data_buffer_size2_[thread_id][depth] = 0;
    }
    return 0;
}

void EmbeddingExplorationEngine::clear_global_queue(int depth)
{
    global_embedding_queue_size_[depth] = 0;
    global_num_ready_embeddings1_[depth] = 0;
    global_num_ready_embeddings2_[depth] = 0;

    global_used_disk_data_buffer_size1_[depth] = 0;
    global_used_disk_data_buffer_size2_[depth] = 0;
}

void EmbeddingExplorationEngine::suspend_thread(const int thread_id, int depth)
{
    {
        std::lock_guard<std::mutex> lk(num_suspended_comp_threads_mutex_[depth]);
        if (global_phases_[depth] != local_phases_[depth][thread_id][0])
        {
            global_phases_[depth] = local_phases_[depth][thread_id][0];
            num_suspended_threads_[depth] = 1;
        }
        else
        {
            num_suspended_threads_[depth] += 1;
        }
        if (num_suspended_threads_[depth] == num_computation_threads_)
        {
            // wake up the main thread
            num_suspended_comp_threads_cv_[depth].notify_one();
        }
    }
    //@@ see https://pubs.opengroup.org/onlinepubs/009696899/functions/pthread_barrier_wait.html
    pthread_barrier_wait(&comp_thread_barrier_[depth]);
    local_phases_[depth][thread_id][0] ^= 1;
}

void EmbeddingExplorationEngine::scatter(id_type id1, id_type id2, int thread_id_, int depth)
{
    int thread_id = thread_id_;
    int socket_id = socket_id_;

    lli required_disk_data_buffer_size1 = (lli)Rrtree_->readNodeLength(id1);
    lli required_disk_data_buffer_size2 = (lli)Srtree_->readNodeLength(id2);

    lli local_queue_size = local_embedding_queue_size_[thread_id][depth + 1];
    lli local_disk_size1 = local_needed_disk_data_buffer_size1_[thread_id][depth + 1];
    lli local_disk_size2 = local_needed_disk_data_buffer_size2_[thread_id][depth + 1];

    if (local_queue_size + 1 > max_num_embeddings_local_queue_ ||
        local_disk_size1 + required_disk_data_buffer_size1 > local_disk_data_buffer_size_ ||
        local_disk_size2 + required_disk_data_buffer_size2 > local_disk_data_buffer_size_)
    {
        // std::cout << "Overflow local mem..." << std::endl;
        int r = flush_local_embeddings(thread_id, depth + 1);
        // @@: r = -1 代表 global_queue 已经满了， 不能再往里面添加新的 embedding 了
        while (r == -1)
        {
            // std::cout << "Overflow global mem..." << std::endl;
            suspend_thread(thread_id, depth);
            // std::cout << "Wakeup... depth = " << depth << std::endl;
            r = flush_local_embeddings(thread_id, depth + 1);
            // std::cout << "flush_local_embeddings done, depth = " << depth << std::endl;
        }
        local_queue_size = local_embedding_queue_size_[thread_id][depth + 1];
        local_disk_size1 = local_needed_disk_data_buffer_size1_[thread_id][depth + 1];
        local_disk_size2 = local_needed_disk_data_buffer_size2_[thread_id][depth + 1];
    }

    // std::cout << "local_queue_size = " << local_queue_size << " " << max_num_embeddings_local_queue_ << std::endl;
    // std::cout << "thread id = " << thread_id << ", next_depth_ = " << next_depth_ << std::endl;
    CompactExtendableEmbedding *e = local_embedding_queues_[thread_id][depth + 1] + local_queue_size;

    // std::cout << "address of e: " << e << std::endl;

    // std::cout << id1 << " " << id2 << std::endl;
    e->id1 = id1;
    e->id2 = id2;
    // * no need to initialize these 2 ptrs
    // e->new_disk_data_1 = (uint8_t *)required_disk_data_buffer_size1;
    // e->new_disk_data_2 = (uint8_t *)required_disk_data_buffer_size2;
    e->disk_data_length_1 = required_disk_data_buffer_size1;
    e->disk_data_length_2 = required_disk_data_buffer_size2;
    // update the queue meta information
    local_embedding_queue_size_[thread_id][depth + 1] = local_queue_size + 1;
    local_needed_disk_data_buffer_size1_[thread_id][depth + 1] = local_disk_size1 + required_disk_data_buffer_size1;
    local_needed_disk_data_buffer_size2_[thread_id][depth + 1] = local_disk_size2 + required_disk_data_buffer_size2;
}

void EmbeddingExplorationEngine::scatter_imbalance(id_type id1,
                                                   id_type id2,
                                                   int thread_id_,
                                                   int depth,
                                                   SimpleRegion const &region)
{
    int thread_id = thread_id_;
    int socket_id = socket_id_;

    lli required_disk_data_buffer_size1 = (lli)(id1 < 0 ? 0 : Rrtree_->readNodeLength(id1));
    lli required_disk_data_buffer_size2 = (lli)(id2 < 0 ? 0 : Srtree_->readNodeLength(id2));

    lli local_queue_size = local_embedding_queue_size_[thread_id][depth + 1];
    lli local_disk_size1 = local_needed_disk_data_buffer_size1_[thread_id][depth + 1];
    lli local_disk_size2 = local_needed_disk_data_buffer_size2_[thread_id][depth + 1];

    if (local_queue_size + 1 > max_num_embeddings_local_queue_ ||
        local_disk_size1 + required_disk_data_buffer_size1 > local_disk_data_buffer_size_ ||
        local_disk_size2 + required_disk_data_buffer_size2 > local_disk_data_buffer_size_)
    {
        // std::cout << "Overflow local mem..." << std::endl;
        int r = flush_local_embeddings(thread_id, depth + 1);
        // @@: r = -1 代表 global_queue 已经满了， 不能再往里面添加新的 embedding 了
        while (r == -1)
        {
            // std::cout << "Overflow global mem..." << std::endl;
            suspend_thread(thread_id, depth);
            // std::cout << "Wakeup... depth = " << depth << std::endl;
            r = flush_local_embeddings(thread_id, depth + 1);
            // std::cout << "flush_local_embeddings done, depth = " << depth << std::endl;
        }
        local_queue_size = local_embedding_queue_size_[thread_id][depth + 1];
        local_disk_size1 = local_needed_disk_data_buffer_size1_[thread_id][depth + 1];
        local_disk_size2 = local_needed_disk_data_buffer_size2_[thread_id][depth + 1];
    }

    CompactExtendableEmbedding *e = local_embedding_queues_[thread_id][depth + 1] + local_queue_size;
    e->id1 = id1;
    e->id2 = id2;
    e->low0 = region.low0;
    e->low1 = region.low1;
    e->high0 = region.high0;
    e->high1 = region.high1;
    e->disk_data_length_1 = required_disk_data_buffer_size1;
    e->disk_data_length_2 = required_disk_data_buffer_size2;
    // update the queue meta information
    local_embedding_queue_size_[thread_id][depth + 1] = local_queue_size + 1;
    local_needed_disk_data_buffer_size1_[thread_id][depth + 1] = local_disk_size1 + required_disk_data_buffer_size1;
    local_needed_disk_data_buffer_size2_[thread_id][depth + 1] = local_disk_size2 + required_disk_data_buffer_size2;
}

void EmbeddingExplorationEngine::scatter_vertex_extendable_embedding(id_type id1, id_type id2)
{
    lli required_disk_data_buffer_size1 = (lli)Rrtree_->readNodeLength(id1);
    lli required_disk_data_buffer_size2 = (lli)Srtree_->readNodeLength(id2);

    // std::cout << "id1 = " << id1 << ", required_disk_data_buffer_size1 = " << required_disk_data_buffer_size1 << std::endl;
    // std::cout << "id2 = " << id2 << ", required_disk_data_buffer_size2 = " << required_disk_data_buffer_size2 << std::endl;

    lli queue_size = global_embedding_queue_size_[0];
    lli used_disk_data_size1 = global_used_disk_data_buffer_size1_[0];
    lli used_disk_data_size2 = global_used_disk_data_buffer_size2_[0];

    if (queue_size + 1 > max_num_embeddings_global_queue_ ||
        used_disk_data_size1 + required_disk_data_buffer_size1 > global_disk_data_buffer_size_ ||
        used_disk_data_size2 + required_disk_data_buffer_size2 > global_disk_data_buffer_size_)
    {
        flush_all_extendable_embeddings();

        queue_size = global_embedding_queue_size_[0];
        used_disk_data_size1 = global_used_disk_data_buffer_size1_[0];
        used_disk_data_size2 = global_used_disk_data_buffer_size2_[0];

        assert(queue_size + 1 <= max_num_embeddings_local_queue_);
        assert(used_disk_data_size1 + required_disk_data_buffer_size1 <= global_disk_data_buffer_size_);
        assert(used_disk_data_size2 + required_disk_data_buffer_size2 <= global_disk_data_buffer_size_);
    }

    CompactExtendableEmbedding *e = global_embedding_queues_[0] + queue_size;
    e->id1 = id1;
    e->id2 = id2;
    // * no need to initialize these 2 ptrs
    // e->new_disk_data_1 = (uint8_t *)required_disk_data_buffer_size1;
    // e->new_disk_data_2 = (uint8_t *)required_disk_data_buffer_size2;
    e->disk_data_length_1 = required_disk_data_buffer_size1;
    e->disk_data_length_2 = required_disk_data_buffer_size2;

    global_embedding_queue_size_[0] = queue_size + 1;
    global_used_disk_data_buffer_size1_[0] += required_disk_data_buffer_size1;
    global_used_disk_data_buffer_size2_[0] += required_disk_data_buffer_size2;
}

void EmbeddingExplorationEngine::spatial_join_func(id_type id1,
                                                   id_type id2,
                                                   NodePtr &n1,
                                                   NodePtr &n2,
                                                   int thread_id,
                                                   int depth)
{
    uint32_t i = 0, j = 0, k;

    if (n1->m_level != n2->m_level)
    {
        std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << std::endl;
        exit(-1);
    }

    // std::cout << "id1 = " << id1 << ", id2 = " << id2 << std::endl;

    // std::cout << "n1->m_children = " << n1->m_children << ", n2->m_children = " << n2->m_children << std::endl;

    while (i < n1->m_children && j < n2->m_children)
    {

        if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0])
        {
            RegionPtr region = n1->m_ptrMBR[i];
            k = j;
            while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0])
            {

                if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] &&
                    region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1])
                {

                    if (n1->m_level == 0 && n2->m_level == 0)
                    {
                        // assert( next_depth_ == 2);
                        // std::cout << "IN SJ....... depth = " << depth << std::endl;
                        count_->add(socket_id_, thread_id, 1);
                        // std::cout << "OUT SJ....... depth = " << depth << std::endl;
                    }
                    else if (n1->m_level == 0)
                    {
                        // std::cout << "aaaaaa next_depth = " << next_depth_ << std::endl;
                        assert(false);
                        scatter(id1, n2->m_pIdentifier[k], thread_id, depth);
                    }
                    else if (n2->m_level == 0)
                    {
                        // std::cout << "bbbbb next_depth = " << next_depth_ << std::endl;
                        // std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << std::endl;
                        assert(false);

                        scatter(n1->m_pIdentifier[i], id2, thread_id, depth);
                    }
                    else
                    {
                        // std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << " " << next_depth_ << std::endl;
                        scatter(n1->m_pIdentifier[i], n2->m_pIdentifier[k], thread_id, depth);
                    }
                }
                k++;
            }
            i++;
        }
        else
        {
            RegionPtr region = n2->m_ptrMBR[j];
            k = i;
            while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0])
            {

                if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] &&
                    region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1])
                {

                    if (n1->m_level == 0 && n2->m_level == 0)
                    {
                        // assert( next_depth_ == 2);
                        // std::cout << "IN SJ....... depth = " << depth << std::endl;
                        count_->add(socket_id_, thread_id, 1);
                        // std::cout << "OUT SJ....... depth = " << depth << std::endl;
                    }
                    else if (n1->m_level == 0)
                    {
                        // std::cout << "aaaaaa next_depth = " << next_depth_ << std::endl;
                        assert(false);
                        scatter(id1, n2->m_pIdentifier[j], thread_id, depth);
                    }
                    else if (n2->m_level == 0)
                    {
                        // std::cout << "bbbbb next_depth = " << next_depth_ << std::endl;
                        // std::cout << "n1->m_level = " << n1->m_level << std::endl;
                        assert(false);
                        scatter(n1->m_pIdentifier[k], id2, thread_id, depth);
                    }
                    else
                    {
                        // std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << " " << next_depth_ << std::endl;
                        scatter(n1->m_pIdentifier[k], n2->m_pIdentifier[j], thread_id, depth);
                    }
                }
                k++;
            }
            j++;
        }
    }
}

void EmbeddingExplorationEngine::spatial_join_func_simplenode_oneside(CompactExtendableEmbedding *e,
                                                                      SimpleNodePtr *n,
                                                                      int thread_id,
                                                                      int depth)
{
    if (e->id1 < 0)
    {
        uint32_t k = 0;
        SimpleRegion region_k = n->returnChildRegion(k);
        while (k < n->m_children && region_k.low0 <= e->high0)
        {
            if (e->low0 <= region_k.high0 &&
                e->low1 <= region_k.high1 &&
                e->high1 >= region_k.low1)
            {
                if (n->m_level == 0)
                {
                    count_->add(socket_id_, thread_id, 1);
                }
                else
                {
                    SimpleRegion ret = region_k.getIntersectingArea(e->low0,
                                                                    e->low1,
                                                                    e->high0,
                                                                    e->high1);
                    scatter_imbalance(-1, n->returnChildIdentifier(k), thread_id, depth, ret);
                }
            }
            k++;
            if (k < n->m_children)
                region_k = n->returnChildRegion(k);
        }
    }
    else if (e->id2 < 0)
    {
        uint32_t k = 0;
        SimpleRegion region_k = n->returnChildRegion(k);
        while (k < n->m_children && region_k.low0 <= e->high0)
        {
            if (e->low0 <= region_k.high0 &&
                e->low1 <= region_k.high1 &&
                e->high1 >= region_k.low1)
            {
                if (n->m_level == 0)
                {
                    count_->add(socket_id_, thread_id, 1);
                }
                else
                {
                    SimpleRegion ret = region_k.getIntersectingArea(e->low0,
                                                                    e->low1,
                                                                    e->high0,
                                                                    e->high1);
                    scatter_imbalance(n->returnChildIdentifier(k), -1, thread_id, depth, ret);
                }
            }
            k++;
            if (k < n->m_children)
                region_k = n->returnChildRegion(k);
        }
    }
    else
        throw std::runtime_error("ERROR...");
}

void EmbeddingExplorationEngine::spatial_join_func_simplenode(id_type id1,
                                                              id_type id2,
                                                              SimpleNodePtr *n1,
                                                              SimpleNodePtr *n2,
                                                              int thread_id,
                                                              int depth)
{
    uint32_t i = 0, j = 0, k;
    while (i < n1->m_children && j < n2->m_children)
    {
        if (n1->returnChildRegion(i).low0 < n2->returnChildRegion(j).low0)
        {
            k = j;
            SimpleRegion region = n1->returnChildRegion(i);
            SimpleRegion region_k = n2->returnChildRegion(k);
            while (k < n2->m_children && region_k.low0 <= region.high0)
            {

                if (region.low1 <= region_k.high1 &&
                    region.high1 >= region_k.low1)
                {

                    if (n1->m_level == 0 && n2->m_level == 0)
                    {
                        count_->add(socket_id_, thread_id, 1);
                    }
                    else if (n1->m_level == 0)
                    {
                        SimpleRegion ret = region.getIntersectingArea(region_k);
                        scatter_imbalance(-1, n2->returnChildIdentifier(k), thread_id, depth, ret);
                    }
                    else if (n2->m_level == 0)
                    {
                        SimpleRegion ret = region.getIntersectingArea(region_k);
                        scatter_imbalance(n1->returnChildIdentifier(i), -1, thread_id, depth, ret);
                    }
                    else
                    {
                        scatter(n1->returnChildIdentifier(i), n2->returnChildIdentifier(k), thread_id, depth);
                    }
                }
                k++;
                if (k < n2->m_children)
                    region_k = n2->returnChildRegion(k);
            }
            i++;
        }
        else
        {
            k = i;
            SimpleRegion region = n2->returnChildRegion(j);
            SimpleRegion region_k = n1->returnChildRegion(k);
            while (k < n1->m_children && region_k.low0 <= region.high0)
            {

                if (region.low1 <= region_k.high1 &&
                    region.high1 >= region_k.low1)
                {

                    if (n1->m_level == 0 && n2->m_level == 0)
                    {
                        count_->add(socket_id_, thread_id, 1);
                    }
                    else if (n1->m_level == 0)
                    {
                        SimpleRegion ret = region.getIntersectingArea(region_k);
                        scatter_imbalance(-1, n2->returnChildIdentifier(j), thread_id, depth, ret);
                    }
                    else if (n2->m_level == 0)
                    {
                        SimpleRegion ret = region.getIntersectingArea(region_k);
                        scatter_imbalance(n1->returnChildIdentifier(k), -1, thread_id, depth, ret);
                    }
                    else
                    {
                        scatter(n1->returnChildIdentifier(k), n2->returnChildIdentifier(j), thread_id, depth);
                    }
                }
                k++;
                if (k < n1->m_children)
                    region_k = n1->returnChildRegion(k);
            }
            j++;
        }
    }
}

// void EmbeddingExplorationEngine::spatial_join_func_simplenode_tmp( id_type id1,
//                                                                 id_type id2,
//                                                                 NodePtr & n1,
//                                                                 NodePtr & n2,
//                                                                 SimpleNodePtr * sn1,
//                                                                 SimpleNodePtr * sn2,
//                                                                 int thread_id,
//                                                                 int depth
//                                                             )
// {
//     uint32_t i = 0, j = 0, k;

//     if (n1->m_level != n2->m_level)
//     {
//         std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << std::endl;
//         exit(-1);
//     }

//     // std::cout << "id1 = " << id1 << ", id2 = " << id2 << std::endl;

//     // std::cout << "n1->m_children = " << n1->m_children << ", n2->m_children = " << n2->m_children << std::endl;

//     while (i < n1->m_children && j < n2->m_children) {

//         if (n1->m_ptrMBR[i]->m_pLow[0] < n2->m_ptrMBR[j]->m_pLow[0]) {
//             RegionPtr region = n1->m_ptrMBR[i];
//             k = j;
//             while (k < n2->m_children && n2->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {

//                 if (region->m_pLow[1] <= n2->m_ptrMBR[k]->m_pHigh[1] &&
//                     region->m_pHigh[1] >= n2->m_ptrMBR[k]->m_pLow[1]) {

//                     if (n1->m_level == 0 && n2->m_level == 0)
//                     {
//                         // assert( next_depth_ == 2);
//                         // std::cout << "IN SJ....... depth = " << depth << std::endl;
//                         count_->add(socket_id_, thread_id, 1);
//                         // std::cout << "OUT SJ....... depth = " << depth << std::endl;
//                     }
//                     else if (n1->m_level == 0)
//                     {
//                         // std::cout << "aaaaaa next_depth = " << next_depth_ << std::endl;
//                         assert(false);
//                         scatter(id1, n2->m_pIdentifier[k], thread_id, depth);
//                     }
//                     else if (n2->m_level == 0)
//                     {
//                         // std::cout << "bbbbb next_depth = " << next_depth_ << std::endl;
//                         // std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << std::endl;
//                         assert(false);

//                         scatter(n1->m_pIdentifier[i], id2, thread_id, depth);
//                     }
//                     else {
//                         // std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << " " << next_depth_ << std::endl;
//                         scatter(n1->m_pIdentifier[i], n2->m_pIdentifier[k], thread_id, depth);
//                     }
//                 }
//                 k++;
//             }
//             i++;
//         }
//         else {
//         //     RegionPtr region = n2->m_ptrMBR[j];
//         //     k = i;
//         //     while (k < n1->m_children && n1->m_ptrMBR[k]->m_pLow[0] <= region->m_pHigh[0]) {

//         //         if (region->m_pLow[1] <= n1->m_ptrMBR[k]->m_pHigh[1] &&
//         //             region->m_pHigh[1] >= n1->m_ptrMBR[k]->m_pLow[1]) {

//         //             if (n1->m_level == 0 && n2->m_level == 0)
//         //             {
//         //                 // assert( next_depth_ == 2);
//         //                 // std::cout << "IN SJ....... depth = " << depth << std::endl;
//         //                 count_->add(socket_id_, thread_id, 1);
//         //                 // std::cout << "OUT SJ....... depth = " << depth << std::endl;
//         //             }
//         //             else if (n1->m_level == 0)
//         //             {
//         //                 // std::cout << "aaaaaa next_depth = " << next_depth_ << std::endl;
//         //                 assert(false);
//         //                 scatter(id1, n2->m_pIdentifier[j], thread_id, depth);
//         //             }
//         //             else if (n2->m_level == 0)
//         //             {
//         //                 // std::cout << "bbbbb next_depth = " << next_depth_ << std::endl;
//         //                 // std::cout << "n1->m_level = " << n1->m_level << std::endl;
//         //                 assert(false);
//         //                 scatter(n1->m_pIdentifier[k], id2, thread_id, depth);
//         //             }
//         //             else {
//         //                 // std::cout << id1 << " " << id2 << " " << n1->m_level << " " << n2->m_level << " " << next_depth_ << std::endl;
//         //                 scatter(n1->m_pIdentifier[k], n2->m_pIdentifier[j], thread_id, depth);
//         //             }
//         //         }
//         //         k++;
//         //     }
//             j++;
//         }
//     }
// }