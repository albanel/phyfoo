cd $1

curl https://www.ensembl.org/info/docs/api/compara/compara_schema.html > generic_file.txt

source ~/beautiful_soup/bin/activate

python read_html.py

deactivate
