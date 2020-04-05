How to get the cusparse BLAS to compile (in Windows)

You need to hack the file cusparse.h (in the folder c:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.2\include\)
Find the pre-processor tags:

#if !defined(_WIN32) (line 7346)

#endif // !defined(_WIN32)  (line 7705)

Out-comment them (in lack of a better way to fix this - nvcc does compile with x64 as default, so I am clueless as to why this win32 thing is an issue)
->

//#if !defined(_WIN32) (line 7346)

//#endif // !defined(_WIN32)  (line 7705)

And then you may proceed with the regular compilation