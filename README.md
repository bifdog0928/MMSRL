# MMSRL
Constructing effective representations of spatial resolved transcriptomics (SRT) data, by appropriately characterizing the coherence in gene expression and histology with the spatial information of each sequencing spot, plays an important role in understanding the organization and function of complex tissues. Although much progress has been made, existing SRT representation methods typically establish local associative relationships for each spot only with those located in its surrounding spatial areas, thus failing to capture longrange correlations between distant regions. In addition, the absence of supervision signals on which cluster (with similar biological functions, pathological states or cell types) each spot should belong to also poses a great challenge in deriving an effective representation of SRT data. To this end, we propose a novel Multi-Modal SRT Representation Learning (MMSRL) method.
![MMSRL Overview](https://github.com/bifdog0928/MMSRL/blob/main/Pic/MMSRL.png?raw=true)

# Requirements

* python~=3.8.19
* numpy~=1.21.5
* numba~=0.55.1
* scanpy~=1.9.3
* torch~=1.11.0
* munkres~=1.1.4
* pandas~=1.3.5
* scikit-learn~=1.0.2
* scikit-misc=~0.2.0
* anndata~=0.8.0
* scipy~=1.7.3
* matplotlib~=3.5.2
* rpy2~=3.4.1
* scvi~=0.6.8

# Datasets

All datasets used in this paper are publicly available and can be downloaded from the links below.
| Datasets | Download link |  
|---------|------------|
| DLPFC | [https://github.com/LieberInstitute/spatialLIBD](https://github.com/LieberInstitute/spatialLIBD) |
| Chicken heart | [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149457) |
| Human breast cancer | [https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0](https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0) | 
| Mouse anterior brain | [https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0) |

# Tutorials

* Perform spatial domain identification on spatial transcriptomic data using MMSRL.
* []{}
