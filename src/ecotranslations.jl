# Getters and setters with microbiome-flavored names

# species -> features
const nfeatures(com::AbstractComMatrix) = nspecies(com)
const getfeature(com::AbstractComMatrix, idx) = getspecies(com, idx)
const featurenames(com::AbstractComMatrix) = specnames(com)
const featuretotals(com::AbstractComMatrix) = speciestotals(com)

# sites -> samples
const nsamples(com::AbstractComMatrix) = nsites(com)
const getsample(com::AbstractComMatrix, idx) = getsite(com, idx)
const samplenames(com::AbstractComMatrix) = sitenames(com)
const sampletotals(com::AbstractComMatrix) = sitetotals(com)
