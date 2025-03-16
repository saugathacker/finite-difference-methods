# **Implementation of EFD, IFD, CNFD**

## **Introduction**

This repository contains an implementation of option pricing estimation using Finite Difference methods in C++ and Python. The code includes methods for solve Explicit Finite Difference Scheme, Implicit Finite Difference Scheme and Crank-Nicolson Finite Difference Scheme.

## **How to Run the Code**

### **C++ Implementation**

Build the project using CMake:

```sh
cmake -Bbuild
cmake --build build
```

Run the compiled executable:

```sh
./build/main
```

### **Python Notebooks**

To run the Jupyter Notebooks:

1. Open a terminal and start Jupyter Notebook:
   ```sh
   jupyter notebook
   ```

## **Dependencies**

### **C++ Requirements:**

- C++20 or later
- CMake
- Standard math libraries (cmath, numbers)

### **Python Requirements:**

- pandas
- numpy
- matplotlib
- scipy
- jupyter

Install required Python packages using:

```sh
pip install pandas numpy matplotlib scipy jupyter
```
