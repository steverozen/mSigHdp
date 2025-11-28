# degenerate_hdpx_extract_components2_many_idential_inputs

    Code
      retvalx
    Output
      $signature
             hdp.1
      1  0.6012711
      2  0.2657635
      3  0.1329654
      4  0.0000000
      5  0.0000000
      6  0.0000000
      7  0.0000000
      8  0.0000000
      9  0.0000000
      10 0.0000000
      
      $signature.post.samp.number
        Signature NumberOfPostSamples
      1     hdp.1                 400
      
      $signature.cdc
         hdp.1
      1   5979
      2    400
      3    399
      4    400
      5    398
      6    400
      7    400
      8    400
      9    400
      10   398
      11   397
      12   398
      13   397
      14   397
      15   396
      16   399
      
      $exposureProbs
            SP.Syn.Abst.1 SP.Syn.Abst.1.1 SP.Syn.Abst.1.2 SP.Syn.Abst.1.3
      hdp.1             1               1               1               1
            SP.Syn.Abst.1.4 SP.Syn.Abst.1.5 SP.Syn.Abst.1.6 SP.Syn.Abst.1.7
      hdp.1               1               1               1               1
            SP.Syn.Abst.1.8 SP.Syn.Abst.1.9 SP.Syn.Abst.1.10 SP.Syn.Abst.1.11
      hdp.1               1               1                1                1
            SP.Syn.Abst.1.12 SP.Syn.Abst.1.13 SP.Syn.Abst.1.14
      hdp.1                1                1                1
      
      $low.confidence.signature
         low confidence hdp.1
      1                     0
      2                     7
      3                     0
      4                     0
      5                     0
      6                     0
      7                     0
      8                     0
      9                     0
      10                    0
      
      $low.confidence.post.samp.number
                   Signature NumberOfPostSamples
      2 low confidence hdp.1                   6
      
      $low.confidence.cdc
         low confidence hdp.1
      1                     7
      2                     0
      3                     0
      4                     0
      5                     0
      6                     0
      7                     0
      8                     0
      9                     0
      10                    0
      11                    3
      12                    1
      13                    1
      14                    2
      15                    0
      16                    0
      
      $extracted.retval
      $extracted.retval$components
            1 2
      1  3595 0
      2  1589 7
      3   795 0
      4     0 0
      5     0 0
      6     0 0
      7     0 0
      8     0 0
      9     0 0
      10    0 0
      
      $extracted.retval$components.post.samples
        Group.1   x
      1       1 400
      2       2   6
      
      $extracted.retval$components.cdc
            1 2
      1  5979 7
      2   400 0
      3   399 0
      4   400 0
      5   398 0
      6   400 0
      7   400 0
      8   400 0
      9   400 0
      10  398 0
      11  397 3
      12  398 1
      13  397 1
      14  397 2
      15  396 0
      16  399 0
      
      $extracted.retval$each.chain.noise.cdc
         summary[[i]]$noise.cdc summary[[i]]$noise.cdc
      1                       3                     11
      2                       0                      0
      3                       0                      1
      4                       0                      0
      5                       0                      2
      6                       0                      0
      7                       0                      0
      8                       0                      0
      9                       0                      0
      10                      0                      2
      11                      0                      0
      12                      0                      1
      13                      0                      2
      14                      1                      0
      15                      1                      3
      16                      1                      0
      
      $extracted.retval$each.chain.noise.clusters
         summary[[i]]$noise.spectrum summary[[i]]$noise.spectrum
      1                            0                           5
      2                            1                           3
      3                            2                           3
      4                            0                           0
      5                            0                           0
      6                            0                           0
      7                            0                           0
      8                            0                           0
      9                            0                           0
      10                           0                           0
      
      $extracted.retval$multi.chains
      Object of class hdpSampleMulti 
       Number of chains: 2 
       Total posterior samples: 400 
       Components: NO. Run hdp_extract_components 
       ----------
       Final hdpState from first chain: 
      Object of class hdpState 
       Number of DP nodes: 16 
       Index of parent DP: 0 1 1 1 1 1 1 1 1 1 ...
       Number of data items per DP: 0 1 1 1 1 1 1 1 1 1 ...
       Index of conparam per DP: 1 2 2 2 2 2 2 2 2 2 ...
       Conparam hyperparameters and current value:
                 Shape Rate      Value
      Conparam 1     1   20 0.07342550
      Conparam 2     1   20 0.05108091
       Number of data categories: 10 
       Number of clusters: 1 
       Initialised with 5 clusters, using random seed 1000044 
      
      $extracted.retval$nsamp
      [1] 400
      
      

