# Parallel out-of-core spatial join processing

## How to Compile

### Dependencies
<ul>
<li> MPI </li> 
<li> NUMA </li>
<li> CMake </li>
</ul>

You might need to manually configure the NUMA library in the CMakeList.txt.
```
cmake -B build
cmake --build build 
```

## Running on generated datasets
### 1. Compile the data generation code
```
g++ gen_data.cpp -o gen_data
```
### 2. Generate 2 datasets to join, each has 1K objects
```
./gen_data 1000 data1.txt
./gen_data 1000 data2.txt
```
### 3. Preprocess 2 datasets into R-trees
```
./build/build_tree data1.txt tree1
./build/build_tree data2.txt tree2
```
### 4. Spatial join tree1 and tree2
```
./build/run_fg tree1 tree2
```
