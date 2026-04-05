# PDB Easter Egg Decorator (PDB Pasqua)

A small [Streamlit](https://streamlit.io/) app that turns [RCSB PDB](https://www.rcsb.org/) structures into pastel Easter egg PNGs from Cα traces (2D PCA projection).

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Run

```bash
streamlit run app.py
```

Enter up to three four-character PDB IDs, click **Hatch my eggs!**, then download each PNG.

## License

Structure data: [RCSB PDB](https://www.rcsb.org/).
