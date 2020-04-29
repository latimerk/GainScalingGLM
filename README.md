This project models Hodgkin-Huxley neurons in response to changes in the stimulus variance.
The simulations are fit with GLMs and the results compared.

Simulations are divided into two parts:

1)gain scaling
    Mease, R. A., Famulare, M., Gjorgjieva, J., Moody, W. J., and Fairhall, A. L. (2013). Emergence of adaptive computation by single neurons in the developing cortex. Journal of Neuroscience

2) Fractional differentiation
    Lundstrom, B. N., Higgs, M. H., Spain, W. J., and Fairhall, A. L. (2008).  Fractional differentiation by neocortical pyramidal neurons. Nature neuroscience


Dependencies: Kenneth's CUDA-based GLM toolbox: https://github.com/latimerk/kGLM
    -Note that GLM fitting is by default run on GPU with CUDA!

        Set path in fitGLMs_LundstromHH.m, fitGLMs_MeaseHH.m
                    measureGLMfits_Lundstrom.m, measureGLMfits_Mease.m

    To compute the Earthmover-distance, we use the following MATLAB function: https://www.mathworks.com/matlabcentral/fileexchange/22962-the-earth-mover-s-distance
        by Ulas Yilmaz

    This function is a wrapper for the EMD computation by Y. Rubner found here http://ai.stanford.edu/~rubner/emd/default.htm  
    and described in Y. Rubner, C. Tomasi, and L. J. Guibas. A Metric for Distributions with Applications to Image Databases. Proceedings of the 1998 IEEE International Conference on Computer Vision, Bombay, India, January 1998, pp. 59-66.

        Set path in quantifyGSdistances.m
