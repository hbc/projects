module purge
source ~/Scripts/loadbipymodules.sh
source  /n/HSPH/local/share/virtualenvs/bipy/bin/activate
cd ~/consults/vv_miRNA-Seq/
python rnaseq_pipeline.py mouse_mapping.full.nohtseqcount.yaml
