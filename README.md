# qp-qiime2
Qiime 2 plugin for Qiita.

Note that this plugin assumes that QIIME 2 is already installed. For instructions follow https://docs.qiime2.org/2017.4/install/.

If you want to add a new plugin, you need to:
* Add the installation plugin in the GitHub Action yml
* Add the new plugin into the `Q2_ALLOWED_PLUGINS` variable in `qp_qiime2.py`
* Add a test to make sure that things work as expected
