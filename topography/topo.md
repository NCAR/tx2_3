# tx2_3v3 Topography Generation

This document describes the workflow used to generate the **tx2_3v3 ocean model topography** and associated files.

---

## Topography Creation Workflow

### 1) Create initial topography

Run:

```
qsub run_create_topo.tx2_3v3.pbs
```

This step creates the basic topography file with **raw cell-level depth statistics** and a **first-cut land mask**.

**Output**

```
topo.sub150.tx2_3v3.SRTM15_V2.4.nc
```

---

### 2) Edit land mask

Notebook:

```
MaskEdit_tx2_3v3.ipynb
```

This notebook edits the land mask, including **removing 1-point bays around Greenland and Antarctica**.

**Output**

```
topo.sub150.tx2_3v3.SRTM15_V2.4.edit1.nc
```

---

### 3) Interpolate and smooth depths

Run:

```
qsub run_append_topo_interp_smooth_tx2_3v3.pbs
```

This step appends a new depth variable **`D_interp`** to the topography file with **lightly smoothed depths**.

---

### 4) Apply manual depth edits

Notebook:

```
Append_topo_edits_tx2_3v3.ipynb
```

This notebook applies **hand edits to depths for straits and channels**.

This step produces the **final model input topography**.

Final topography variable used by the model:

```
D_edit2
```

---

### 5) Create channel width file

Notebook:

```
Channel_width_tx2_3v3.ipynb
```

Creates the **channel width file** used by the model.

---

## Optional Steps (not used in v3)

### 6) Channel depth edits for runtime experiments

Notebook:

```
Channel_topo_tx2_3v3.ipynb
```

This notebook allows **local topography depth edits** for runtime experiments.

---

### 7) Incorporate channel edits into the base topography

If changes from step 6 must be incorporated into the topography generated in step 3, run:

```
Append_topo_edits.ipynb
```

---

## Validation / Diagnostic Notebooks

Several notebooks in the source directory are used to verify each step of the workflow.

### CompareFirstMasks.ipynb

Ensures **step 1 reproduces the tx2_3v2 mask**.

### CompareEditedMasks.ipynb

Checks the mask after step 2 to confirm that **only intended edits were made**.

### CheckRunoffPoints.ipynb

Similar validation focused on **locations where runoff points were intentionally modified**.

### CompareDepths.ipynb

Compares the **final depths from step 4** with the **input and output files from tx2_3v2**.

Expected behavior:

* Depth differences should appear **only at the ~30 points modified for tx2_3v3**.
* Small roundoff differences (**O(10⁻³ m)**) relative to tx2_3v2 may occur, possibly due to the **new grid file**.

