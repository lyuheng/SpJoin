#pragma once
//algo: implemented according to http://www.parallellabs.com/2010/10/25/practical-concurrent-queue-algorithm

#include <mutex>
#include <atomic>
#include <condition_variable>

#include <iostream>

template <typename T> 
struct node_t
{
    T value;
    std::atomic<node_t<T> *> next;
};

template <typename T> 
class conque
{
public:
    typedef node_t<T> NODE;
    NODE *head; //locked by q_h_lock
    NODE *tail;
    std::mutex q_h_lock;
    std::mutex q_t_lock;
    
    conque() {
        head = tail = new NODE;
        head->next = NULL;
    }
    
    conque(const conque<T>& q) {
    	//required due to ReqQueue's initialization calls resize(.)
    	//the trouble is caused by mutex-elements, which cannot be copied
		//the function is not actually called, we need to define it just because resize(.) may copy
    	std::cout << "conque's copy constructor called, should never happen !!!" << std::endl;
    	exit(-1);
		head = tail = new NODE;
		head->next = NULL;
	}

    ~conque() // prevent memory leak 
    {
        delete head;
    }

    void enqueue(T value) {
        NODE *node = new NODE;
        node->value = value;
        node->next = NULL;
        std::unique_lock<std::mutex> lck(q_t_lock);
        tail->next.store(node, std::memory_order_relaxed); //atomic "store"
        tail = node;
    }
    
    bool dequeue(T & to_get) {
        std::unique_lock<std::mutex> lck(q_h_lock);
        NODE *node = head;
        NODE *new_head = node->next.load(std::memory_order_relaxed); //atomic "load"
        if(new_head == NULL) //queue is empty
        {
        	lck.unlock();
        	return false;
        }
        to_get = new_head->value;
        head = new_head;

        lck.unlock();
        delete node;
        return true;
    }

    bool empty()
    {
    	std::unique_lock<std::mutex> lck(q_h_lock);
		NODE *node = head->next.load(std::memory_order_relaxed); //atomic "load"
		lck.unlock();
		if(node == NULL) return true;
		else return false;
    }
};