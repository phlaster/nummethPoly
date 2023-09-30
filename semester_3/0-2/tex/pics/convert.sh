for file in *.svg; do
    inkscape -D "$file" -o "${file%.svg}.pdf" --export-latex
done
