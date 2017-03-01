# Create symlink if necessary
if [ ! -f dropbox ]; then
    ln -s ~/Dropbox/Shared/HBC/sahay_klf9_overexpression_rnaseq dropbox
fi

cp -r data-raw dropbox
cp -r results dropbox
cp -r *_files dropbox

cp *.bib dropbox
cp *.html dropbox
cp *.pdf dropbox
cp *.R dropbox
cp *.Rmd dropbox
cp *.yaml dropbox
