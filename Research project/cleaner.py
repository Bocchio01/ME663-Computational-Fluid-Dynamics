import os

def delete_tex_aux_files(folder):
    extensions = [
        '.aux', '.log', '.out', '.toc', '.bbl', '.blg', '.synctex.gz', '.fdb_latexmk',
        '.fls', '.gz', '.idx', '.ilg', '.ind', '.lof', '.lot', '.maf', '.mtc',
        '.mtc0', '.nav', '.nlo', '.nls', '.snm', '.spl', '.synctex.gz', '.tdo', '.vrb',
        '.xdy', '.acn', '.acr', '.alg', '.glg', '.glo', '.gls', '.ist', '.loa', '.lof',
        '.bcf', '.run.xml', '.synctex(busy)', '.synctex', '.synctex.gz(busy)', '.synctex.gz', '.synctex.gz(busy)',
    ]

    for root, _, files in os.walk(folder):
        for file in files:
            if any(file.endswith(ext) for ext in extensions):
                file_path = os.path.join(root, file)
                os.remove(file_path)
                print(f"Deleted: {file_path}")

if __name__ == "__main__":
    folder_to_clean = "./"

    delete_tex_aux_files(folder_to_clean)
