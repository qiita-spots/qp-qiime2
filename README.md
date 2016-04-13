# qp-qiime2
Qiime 2 plugin for Qiita

To install and use this see the following steps. As this is depending on several experimental packages, it's possible that these steps won't work, but we'll make an effort to keep them up to date.

```bash
conda create -n qp-qiime2 python=3.5 scikit-bio
source activate qp-qiime2
pip install https://github.com/biocore/qiime2/archive/master.zip https://github.com/qiime2/q2-feature-table/archive/master.zip https://github.com/gregcaporaso/qp-qiime2/archive/master.zip
# the a b c here are meaningless right now, just requirements of the interface
start-qiime2 a b c
```
