# Getters and setters with microbiome-flavored names

# species -> features
nfeatures(com::AbstractComMatrix) = nspecies(com)
getfeature(com::AbstractComMatrix, idx) = getspecies(com, idx)
featurenames(com::AbstractComMatrix) = specnames(com)
featuretotals(com::AbstractComMatrix) = speciestotals(com)

# sites -> samples
nsamples(com::AbstractComMatrix) = nsites(com)
getsample(com::AbstractComMatrix, idx) = getsite(com, idx)
samplenames(com::AbstractComMatrix) = sitenames(com)
sampletotals(com::AbstractComMatrix) = sitetotals(com)
