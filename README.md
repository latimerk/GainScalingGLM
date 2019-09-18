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

    Earthmover-distance: https://www.mathworks.com/matlabcentral/fileexchange/22962-the-earth-mover-s-distance
        by Ulas Yilmaz

        Set path in quantifyGSdistances.m
