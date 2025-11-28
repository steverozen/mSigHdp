# degenerate_hdpx_extract_components

    Code
      retvalx
    Output
      $signature
      numeric(0)
      
      $signature.post.samp.number
      [1] Group.1 x      
      <0 rows> (or 0-length row.names)
      
      $signature.cdc
       
      1
      2
      3
      
      $exposureProbs
           [,1]
      
      $low.confidence.signature
         low confidence hdp.1
      1                   800
      2                     0
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
      1 low confidence hdp.1                   1
      
      $low.confidence.cdc
        low confidence hdp.1
      1                  800
      2                  400
      3                  400
      
      $extracted.retval
      $extracted.retval$components
           1
      1  800
      2    0
      3    0
      4    0
      5    0
      6    0
      7    0
      8    0
      9    0
      10   0
      
      $extracted.retval$components.post.samples
        Group.1 x
      1       1 1
      
      $extracted.retval$components.cdc
          1
      1 800
      2 400
      3 400
      
      $extracted.retval$each.chain.noise.cdc
      data frame with 0 columns and 3 rows
      
      $extracted.retval$each.chain.noise.clusters
      data frame with 0 columns and 10 rows
      
      $extracted.retval$multi.chains
      Object of class hdpSampleMulti 
       Number of chains: 2 
       Total posterior samples: 400 
       Components: NO. Run hdp_extract_components 
       ----------
       Final hdpState from first chain: 
      Object of class hdpState 
       Number of DP nodes: 3 
       Index of parent DP: 0 1 1 ...
       Number of data items per DP: 0 1 1 ...
       Index of conparam per DP: 1 2 2 ...
       Conparam hyperparameters and current value:
                 Shape Rate       Value
      Conparam 1     1   20 0.130470280
      Conparam 2     1   20 0.005143543
       Number of data categories: 10 
       Number of clusters: 1 
       Initialised with 5 clusters, using random seed 1000044 
      
      $extracted.retval$nsamp
      [1] 400
      
      

# degenerate_hdpx_extract_components_warnings

    Unable to estimate exposures; no exposure-related output generated

