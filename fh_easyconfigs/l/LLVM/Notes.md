# LLVM build notes

#### limits.h missing

Compile error
```
/build/LLVM/11.1.0/GCCcore-11.2.0/llvm-11.1.0.src/utils/benchmark/src/benchmark_register.h:17:30: error: numeric_limits is not a member of std
   17 |   static const T kmax = std::numeric_limits<T>::max();
      |                              ^~~~~~~~~~~~~~
```

#### Fix
```
include <limits>
```
LLVM/11.1.0/GCCcore-11.2.0/llvm-11.1.0.src/utils/benchmark/src/benchmark_register.h:17:30: error: numeric_limits is not a member of std

LLVM-11.1.0-GCCcore-11.2.0.eb
