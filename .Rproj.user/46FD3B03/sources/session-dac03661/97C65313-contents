import os
import shutil

# Target folders
folders = {
    "scripts": [".R", ".Rproj", ".py", ".ipynb", ".sh"],
    "docs": [".docx", ".pdf", ".md", ".txt"],
    "results": [".png", ".jpg", ".jpeg", ".pdf", ".svg", ".xlsx"],
    "data": [
        ".vcf", ".vcf.gz", ".tbi", ".csi", ".bcf",
        ".bed", ".bim", ".fam", ".frq", ".raw", ".clst",
        ".eigenvec", ".eigenval", ".panel", ".range", ".index"
    ]
}

# Ensure folders exist
for folder in folders:
    os.makedirs(folder, exist_ok=True)

# Move files
for file in os.listdir("."):
    if os.path.isfile(file):  
        ext = os.path.splitext(file)[1]  # keep case sensitivity
        moved = False

        # Special handling for compressed files
        if file.endswith(".vcf.gz") or file.endswith(".vcf.gz.tbi") or file.endswith(".vcf.gz.csi"):
            shutil.move(file, os.path.join("data", file))
            print(f"Moved {file} → data/")
            moved = True
        else:
            for folder, extensions in folders.items():
                if ext in extensions:
                    shutil.move(file, os.path.join(folder, file))
                    print(f"Moved {file} → {folder}/")
                    moved = True
                    break

        if not moved:
            print(f"Skipped {file} (no matching rule)")
