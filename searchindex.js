Search.setIndex({"docnames": ["README", "initial_conditions/ics", "maximum_thickness/max_h", "mesh/mesh", "runoff_mapping/runoff", "salinity_restoring/sss_restoring", "supergrid/orca_gridgen", "tidal_energy_dissipation/regrid_tidal_dissipation", "topography/Append_topo_edits", "topography/Channel_topo_tx2_3v2", "topography/Channel_width_tx2_3v2", "topography/MaskEdit_tx2_3v2"], "filenames": ["README.md", "initial_conditions/ics.md", "maximum_thickness/max_h.ipynb", "mesh/mesh.md", "runoff_mapping/runoff.md", "salinity_restoring/sss_restoring.md", "supergrid/orca_gridgen.md", "tidal_energy_dissipation/regrid_tidal_dissipation.ipynb", "topography/Append_topo_edits.ipynb", "topography/Channel_topo_tx2_3v2.ipynb", "topography/Channel_width_tx2_3v2.ipynb", "topography/MaskEdit_tx2_3v2.ipynb"], "titles": ["2/3 degree global configuration", "Introduction", "Define maximum layer thicknesses H_max(z)", "Generate an ESMF mesh", "Runoff mapping files", "Introduction", "Supergrid", "Wave Dissipation dataset", "Append hand edit topography modifications to the interpolated topography", "Hand edit topography on selected sills and straits", "Set sub-grid widths for selected channels", "Ocean Mask Edits for tx2_3v2"], "terms": {"thi": [0, 1, 2, 4, 5, 6, 7, 11], "repositori": 0, "document": 0, "all": [0, 4, 11], "step": [0, 2, 4], "need": [0, 6, 7], "creat": [0, 2, 3, 6, 7, 10, 11], "nomin": [0, 3, 6, 8, 9, 10, 11], "us": [0, 2, 3, 6, 8, 9], "cesm": [0, 2, 4], "mom6": [0, 2, 6, 7, 8, 9, 10, 11], "first": [0, 7], "clone": 0, "your": [0, 1, 2, 3, 5], "account": 0, "own": 0, "branch": 0, "work": [0, 4, 7], "git": 0, "checkout": 0, "b": [0, 2], "nameofyourbranch": 0, "conda": [0, 7, 8, 9], "env": [0, 7, 8, 9], "f": [0, 10], "yml": 0, "activ": 0, "environment_nam": 0, "an": [1, 4, 5, 6], "overview": [1, 5], "project": [1, 5], "feel": [1, 5], "free": [1, 5], "updat": [1, 2, 5, 8], "accordingli": [1, 2, 4, 5], "from": [2, 4, 7, 8, 9, 10, 11], "mom6_tool": 2, "m6plot": 2, "import": [2, 7, 8, 9, 10, 11], "xyplot": 2, "mystat": 2, "xarrai": [2, 7, 8, 9, 10, 11], "xr": [2, 7, 8, 9, 10, 11], "numpi": [2, 7, 8, 9, 10, 11], "np": [2, 7, 8, 9, 10, 11], "matplotlib": [2, 7, 8, 9, 10, 11], "rcparam": 2, "font": 2, "size": [2, 7, 9], "16": [2, 8, 9, 10, 11], "color": [2, 7, 10], "pyplot": [2, 8, 9, 10, 11], "plt": [2, 8, 9, 10, 11], "scipi": 2, "optim": 2, "curve_fit": 2, "datetim": [2, 7, 8, 9, 10, 11], "inlin": [2, 8, 9, 10, 11], "The": [2, 4, 6, 9], "follow": [2, 4, 6, 11], "paramet": [2, 4], "must": [2, 4, 6], "set": [2, 3, 4], "case": [2, 8, 9, 10], "name": [2, 4, 8, 9, 10, 11], "chang": [2, 6], "each": [2, 11], "configur": [2, 6, 8, 9], "grid_nam": 2, "tx2_3v2": [2, 8, 9, 10], "path": [2, 8, 9, 10, 11], "grid": [2, 3, 4, 6, 8, 9, 11], "baselin": 2, "add": [2, 4], "email": 2, "address": 2, "author": [2, 7], "gustavo": [2, 8], "marqu": [2, 8], "gmarqu": [2, 4, 8], "ucar": [2, 7, 8, 9, 10, 11], "edu": [2, 7, 8, 9, 10, 11], "max_layer_thick": 2, "gfdl": 2, "s": [2, 7, 10], "om4": 2, "comparison": [2, 6], "purpos": 2, "max_h": 2, "arrai": [2, 8, 9], "400": [2, 8], "0": [2, 4, 7, 8, 9, 10, 11], "409": [2, 9, 11], "63": [2, 11], "410": [2, 9], "32": [2, 9, 10, 11], "75": [2, 9, 10, 11], "411": [2, 9], "07": 2, "52": [2, 11], "7": [2, 7, 8, 9, 10, 11], "86": [2, 8, 9, 10, 11], "412": [2, 9], "13": [2, 9, 10, 11], "24": [2, 9], "35": [2, 9, 10], "45": [2, 8, 11], "54": [2, 11], "71": [2, 8, 9, 10, 11], "79": [2, 8, 11], "93": [2, 11], "413": [2, 9, 11], "06": 2, "12": [2, 8, 9, 10, 11], "18": [2, 8, 9, 11], "29": [2, 8, 9, 10, 11], "34": [2, 11], "39": [2, 8, 9, 10, 11], "44": [2, 9, 10, 11], "49": [2, 8, 11], "58": [2, 11], "62": [2, 9, 11], "67": [2, 8, 9, 10, 11], "78": [2, 8, 9, 11], "82": [2, 8, 9, 11], "9": [2, 9, 10, 11], "97": [2, 8], "414": [2, 9], "03": 2, "1": [2, 4, 6, 7, 8, 9, 10, 11], "19": [2, 6, 8, 9, 11], "22": [2, 9, 10, 11], "27": [2, 11], "3": [2, 4, 6, 7, 8, 9, 10, 11], "33": [2, 8, 9, 10, 11], "38": [2, 9, 10, 11], "41": [2, 8, 9, 10, 11], "43": [2, 10, 11], "46": [2, 8, 9, 10, 11], "48": [2, 11], "51": [2, 8, 9, 10, 11], "53": [2, 7, 11], "55": [2, 11], "6": [2, 8, 9, 10, 11], "65": [2, 9, 11], "69": [2, 8, 9, 10, 11], "73": [2, 8, 9, 10, 11], "77": [2, 9, 11], "83": [2, 11], "glade": [2, 4, 7, 8, 9, 10, 11], "derecho": 2, "scratch": [2, 4], "g": 2, "e23_a16g": 2, "gjrav4": 2, "tl319_t232_hycom1_n75": 2, "2024": [2, 8, 9, 10, 11], "hycom1_exploration_n75": 2, "run": [2, 3, 4, 11], "grd": 2, "open_dataset": [2, 7, 8, 9, 10, 11], "h": [2, 7], "static": 2, "nc": [2, 3, 4, 6, 7, 8, 9, 10, 11], "load": [2, 3, 6], "interfac": 2, "zstar": 2, "dz": 2, "vertical_coordin": 2, "eta": 2, "lt": [2, 8, 9, 10, 11], "dataarrai": 2, "x27": [2, 8, 9, 10, 11], "gt": [2, 8, 9, 10, 11], "250000e": 2, "00": 2, "750000e": 2, "8": [2, 6, 7, 8, 9, 10, 11], "125000e": 2, "01": 2, "375000e": 2, "625000e": 2, "875000e": 2, "876372e": 2, "130488e": 2, "387348e": 2, "646951e": 2, "909299e": 2, "017439e": 2, "02": [2, 8, 9, 10, 11], "044253e": 2, "071894e": 2, "101483e": 2, "134005e": 2, "170264e": 2, "211252e": 2, "258301e": 2, "313104e": 2, "377775e": 2, "455000e": 2, "548061e": 2, "661100e": 2, "799286e": 2, "968920e": 2, "177837e": 2, "435695e": 2, "754256e": 2, "147966e": 2, "634415e": 2, "234923e": 2, "975393e": 2, "887171e": 2, "007975e": 2, "383154e": 2, "006720e": 2, "212538e": 2, "463554e": 2, "769054e": 2, "140079e": 2, "589734e": 2, "133543e": 2, "789850e": 2, "580300e": 2, "505995e": 2, "coordin": [2, 6, 8, 9, 10, 11], "float64": [2, 8, 9, 10, 11], "25": [2, 8, 9, 10, 11], "79e": 2, "58e": 2, "506e": 2, "attribut": [2, 8, 9, 10, 11], "long_nam": [2, 8, 9], "rho": 2, "unit": [2, 8, 9], "meter": [2, 4, 8, 9, 10, 11], "cartesian_axi": 2, "posit": 2, "upxarrai": 2, "751": 2, "11": [2, 8, 9, 10, 11], "134e": 2, "03arrai": 2, "float641": 2, "03long_nam": 2, "rhounit": 2, "metercartesian_axi": 2, "zposit": 2, "uparrai": 2, "up": 2, "hybrid": 2, "target": 2, "densiti": 2, "vgrid": 2, "hybrid_75layer_zstar2": 2, "50m": 2, "2020": 2, "23": [2, 9, 11], "sigma2": 2, "input": [2, 4, 7], "fname": 2, "mom_ic": 2, "drop": 2, "time": [2, 11], "hist_nam": 2, "nativ": 2, "0001": 2, "ds": 2, "isel": [2, 7], "def": 2, "get_quantil": 2, "da": 2, "95": [2, 11], "dim": 2, "yh": 2, "xh": 2, "threshold": 2, "filter": 2, "find": 2, "greater": 2, "than": 2, "equal": 2, "da_quantil": 2, "where": [2, 4, 7, 11], "fals": [2, 4], "top_5_percent_valu": 2, "dropna": 2, "how": 2, "ani": [2, 6], "return": 2, "fig": [2, 11], "ax": [2, 8, 9, 10, 11], "subplot": [2, 11], "figsiz": [2, 11], "10": [2, 8, 9, 10, 11], "k": [2, 3, 7, 10], "rang": [2, 8, 9], "len": 2, "data": 2, "flatten": 2, "els": [2, 6], "sum": [2, 11], "skipna": 2, "sc": 2, "scatter": 2, "ones": [2, 8], "c": 2, "norm": [2, 7], "lognorm": [2, 7], "vmin": [2, 7, 11], "vmax": [2, 7, 11], "6000": 2, "cmap": 2, "rdylbu_r": 2, "plot": [2, 7, 10, 11], "label": [2, 11], "dz_max": 2, "legend": 2, "invert_yaxi": 2, "set_xlim": [2, 11], "450": 2, "colorbar": 2, "xlabel": 2, "ylabel": 2, "index": [2, 8, 9, 10, 11], "lath": [2, 8, 9, 10, 11], "lonh": [2, 8, 9, 10, 11], "titl": [2, 4, 8, 9], "90th": 2, "bin": 2, "arang": [2, 8], "200": [2, 11], "mean_90": 2, "zero": 2, "std_90": 2, "20": [2, 11], "data1": 2, "ma": 2, "masked_invalid": 2, "wet": 2, "min": 2, "area_t": 2, "nrow": 2, "ncol": [2, 11], "title1": 2, "format": [2, 7, 8, 10, 11], "masked_wher": 2, "geolon": [2, 8, 9, 10, 11], "geolat": [2, 8, 9, 10, 11], "clim": 2, "50": [2, 9, 11], "axi": 2, "nbin": 2, "sigma": 2, "title2": 2, "data2": 2, "subplots_adjust": 2, "top": 2, "hist": 2, "title3": 2, "histogram": 2, "set_titl": 2, "tmp": 2, "ipykernel_2588": 2, "956470891": 2, "py": [2, 3, 7, 8, 9], "runtimewarn": 2, "more": 2, "figur": 2, "have": 2, "been": [2, 6, 8, 9], "open": [2, 10], "through": 2, "ar": [2, 11], "retain": 2, "until": 2, "explicitli": 2, "close": 2, "mai": 2, "consum": 2, "too": 2, "much": 2, "memori": 2, "To": [2, 6], "control": 2, "warn": [2, 7, 8, 9, 11], "see": [2, 3, 4], "max_open_warn": 2, "49999987": 2, "64398286": 2, "52906127": 2, "528046": 2, "17002452": 2, "82966101": 2, "85044182": 2, "7335062": 2, "78472144": 2, "91840888": 2, "14788992": 2, "91947574": 2, "95183546": 2, "549369": 2, "73346407": 2, "57874192": 2, "15": [2, 9, 10, 11], "82926689": 2, "14": [2, 9, 10, 11], "50177195": 2, "41650331": 2, "4654436": 2, "72917385": 2, "09291042": 2, "78575757": 2, "17": [2, 9, 11], "89876773": 2, "21": [2, 9], "03573048": 2, "23213729": 2, "30": [2, 8, 9, 10, 11], "03095485": 2, "55386447": 2, "90016546": 2, "36": [2, 8, 9, 10, 11], "00712351": 2, "37": [2, 9], "10292082": 2, "87783535": 2, "56": [2, 8, 9, 10, 11], "20806657": 2, "58357165": 2, "72": [2, 8, 9, 10, 11], "15017817": 2, "80": [2, 8, 9, 10, 11], "16838769": 2, "94": 2, "40601621": 2, "116": [2, 9, 10, 11], "13228634": 2, "114": 2, "70040409": 2, "129": [2, 11], "6036234": 2, "131": 2, "76797368": 2, "2640204": 2, "109": 2, "06462448": 2, "111": [2, 11], "97581924": 2, "117": [2, 11], "84792669": 2, "124": 2, "64797682": 2, "139": 2, "74287665": 2, "157": 2, "62992308": 2, "182": 2, "63635534": 2, "218": 2, "98859406": 2, "232": 2, "00104179": 2, "215": 2, "30296714": 2, "201": 2, "19663676": 2, "196": [2, 11], "13597868": 2, "212": 2, "3738024": 2, "213": 2, "79715404": 2, "245": [2, 9, 10, 11], "87368957": 2, "99288766": 2, "247": 2, "73433881": 2, "264": [2, 8], "28649986": 2, "281": [2, 8, 9, 10, 11], "17038354": 2, "300": [2, 8, 11], "34660854": 2, "345": [2, 11], "52134604": 2, "329": [2, 11], "40704216": 2, "356": 2, "7702918": 2, "386": [2, 9], "4174176": 2, "389": [2, 9], "15544909": 2, "392": [2, 9], "40018901": 2, "05565569": 2, "135": 2, "35851589": 2, "66796584": 2, "0637156": 2, "01814795": 2, "03510754": 2, "tanh_funct": 2, "x": [2, 4], "d": [2, 7, 9, 10, 11], "note": [2, 3, 4, 7], "extra": 2, "y": [2, 7, 9, 10, 11], "popt": 2, "pcov": 2, "extract": 2, "gener": [2, 4, 6, 8, 9, 10, 11], "y_fit": 2, "origin": [2, 6, 8, 9], "r": [2, 7, 8, 9, 10], "show": 2, "print": [2, 8, 9, 10, 11], "187": 2, "08972100760633": 2, "056399678180667163": 2, "80764670977144": 2, "189": 2, "62178708778868": 2, "dz_max_90": 2, "88978622": 2, "05124971": 2, "2318273": 2, "43375969": 2, "65954471": 2, "91196546": 2, "19412109": 2, "50946059": 2, "86181952": 2, "25546": 2, "69511391": 2, "18602944": 2, "73402094": 2, "34552182": 2, "02764048": 2, "78821859": 2, "63589128": 2, "58014822": 2, "63139452": 2, "80100975": 2, "10140314": 2, "54606232": 2, "14959251": 2, "92774221": 2, "89741083": 2, "07663286": 2, "26": [2, 8, 9, 10, 11], "4845323": 2, "14124041": 2, "06776905": 2, "28583156": 2, "81760291": 2, "42": [2, 8, 9, 10, 11], "68541118": 2, "91135357": 2, "51683185": 2, "52200507": 2, "61": [2, 8, 9, 10, 11], "94516131": 2, "80201558": 2, "74": [2, 9, 11], "1049473": 2, "86219824": 2, "88": [2, 8, 9, 10, 11], "07705989": 2, "74708725": 2, "103": [2, 11], "86338297": 2, "112": 2, "41000134": 2, "121": [2, 11], "36352345": 2, "130": [2, 11], "69285331": 2, "140": 2, "35927768": 2, "150": 2, "31682062": 2, "160": 2, "51290687": 2, "170": [2, 9, 11], "8893278": 2, "181": 2, "38348107": 2, "191": [2, 11], "92983316": 2, "202": 2, "46153511": 2, "91210815": 2, "223": [2, 11], "21710955": 2, "233": [2, 11], "31569094": 2, "243": [2, 11], "15197092": 2, "252": 2, "67615996": 2, "261": [2, 8], "84539616": 2, "270": [2, 8, 11], "62427305": 2, "278": [2, 8], "98506237": 2, "286": [2, 8, 9, 10, 11], "90765387": 2, "294": [2, 8, 11], "37924899": 2, "301": [2, 9], "39385482": 2, "307": 2, "95162933": 2, "314": [2, 9, 11], "05812899": 2, "319": [2, 9], "72350586": 2, "324": [2, 9, 11], "96169559": 2, "78962955": 2, "334": [2, 9, 11], "22649664": 2, "338": [2, 11], "29307219": 2, "342": [2, 9, 11], "01112468": 2, "40290465": 2, "348": [2, 11], "49071585": 2, "351": [2, 11], "29656501": 2, "353": [2, 8, 9, 11], "84188435": 2, "dataset": [2, 10, 11], "coord": [2, 7], "attr": [2, 7, 8, 9, 11], "018": 2, "now": [2, 7, 8], "strftime": [2, 7, 9, 10, 11], "date_fmt": 2, "dz_max_": 2, "_": [2, 3, 8, 9, 10], "to_netcdf": [2, 7, 8, 9, 11], "supergrid": 3, "ocean": [3, 4, 6, 8, 9, 10], "mask": [3, 4, 8, 9, 10], "topographi": [3, 11], "netcdf": [3, 4, 6, 7], "file": [3, 6], "model": [3, 6, 8, 9, 10, 11], "python": 3, "gen_nc_grid": 3, "make": [3, 4, 6, 7], "sure": [3, 4, 7], "variabl": [3, 4, 8, 9, 10, 11], "nc_topo": 3, "point": [3, 4, 6, 8, 9, 10, 11], "right": 3, "convert": 3, "scrip": [3, 4], "convent": 3, "modul": [3, 6], "ncl": 3, "gen_scrip": 3, "ncgrdfilepath": 3, "maskfilepath": 3, "dstgridpath": 3, "qsub": 3, "create__mesh": 3, "file_i": 3, "file_o": 3, "If": 3, "script": 3, "successfulli": 3, "you": 3, "should": [3, 6], "tx": 3, "_mesh_yymmdd": 3, "directori": 3, "cdf5": [3, 7], "nccopi": [3, 7], "tx2_3_mesh_yymmdd": 3, "tx2_3_mesh_yymmdd_cdf5": 3, "In": 4, "section": 4, "we": [4, 6], "generag": 4, "spread": 4, "liquid": 4, "frozen": 4, "water": [4, 11], "cime": 4, "provid": 4, "tool": 4, "help": 4, "pleas": [4, 6, 8, 9], "build": 4, "runoff_map": 4, "execut": [4, 6], "instruct": 4, "here": 4, "befor": 4, "proceed": 4, "next": [4, 6], "exampl": 4, "namelist": 4, "tx2_3": [4, 6, 8], "map_jra_to_tx2_3": 4, "nml": 4, "input_nml": 4, "gridtyp": 4, "ob": [4, 8, 9, 10, 11], "file_roff": 4, "p": 4, "omwg_dev": 4, "jra55": 4, "domain": 4, "190212": 4, "roff": 4, "jra025m": 4, "190213": 4, "file_ocn": 4, "mesh": [4, 7], "tx2_3_scrip_yymmdd": 4, "file_nn": 4, "map_jra_to_tx2_3_nn": 4, "yymmdd": 4, "file_smooth": 4, "map_jra_to_tx2_3_sm_e333r100_yymmdd": 4, "file_new": 4, "map_jra_to_tx2_3_nnsm_e333r100_yymmdd": 4, "nearest": 4, "neighbor": 4, "smooth": 4, "efold": 4, "1000000": 4, "rmax": 4, "300000": 4, "step1": 4, "true": [4, 7, 8, 9, 11], "step2": 4, "step3": 4, "modifi": [4, 6], "can": [4, 6], "divid": 4, "four": 4, "categori": 4, "type": [4, 6], "rtm": 4, "720": [4, 7], "360": [4, 11], "ascii": 4, "xc": 4, "yc": 4, "xv": 4, "yv": 4, "area": [4, 8, 9, 10, 11], "contain": 4, "grid_area": 4, "along": 4, "typic": 4, "rdirc": 4, "OR": 4, "cell": [4, 8, 9, 10, 11], "below": 4, "file_ocn_coastal_mask": 4, "onli": [4, 6], "coastal": 4, "specifi": 4, "file_ocn_coast_mask": 4, "standard": 4, "includ": [4, 6], "distanc": 4, "default": [4, 6], "maximum": [4, 8, 9, 10, 11], "radiu": 4, "effect": 4, "string": 4, "unset": 4, "restrict_smooth_src_to_nn_dest": 4, "option": 4, "limit": 4, "sourc": [4, 8, 9, 10, 11], "just": 4, "get": 4, "instead": [4, 8, 9], "comput": 4, "multipli": 4, "two": 4, "togeth": 4, "output": 4, "field": [4, 7], "nn": 4, "smoother": 4, "combin": 4, "nnsm": 4, "jra_to_ocn": 4, "sh": 4, "sandbox": 4, "cesm2_3_beta08_carib": 4, "gen_mapping_fil": 4, "runoff_to_ocn": 4, "src": [4, 8, 9, 10, 11], "env_mach_specif": 4, "export": 4, "tmpdir": 4, "map_r05_to_tx2_3": 4, "cseg": 4, "inputdata": 4, "lnd": 4, "clm2": 4, "rtmdata": 4, "05": 4, "061026": 4, "tx2_3_scrip_230415": 4, "map_r05_to_tx2_3_nn": 4, "230415": 4, "map_r05_to_tx2_3_sm_e1000r1000_230415": 4, "map_r05_to_tx2_3_nnsm_e1000r1000_230415": 4, "r05_to_ocn": 4, "gland4km": 4, "map_greenland_4km_to_tx2_3": 4, "glc": 4, "cism": 4, "griddata": 4, "scripgrid_greenland_4km_epsg3413_c170414": 4, "map_gland4km_to_tx2_3_nn": 4, "map_gland4km_to_tx2_3_sm_e1000r300_230415": 4, "map_gland4km_to_tx2_3_nnsm_e1000r300_230415": 4, "gland4km_to_ocn": 4, "orca_gridgen": 6, "which": 6, "reli": 6, "nemo": 6, "framework": 6, "tripolar": 6, "For": [6, 11], "complet": 6, "descript": [6, 7, 8, 9, 10, 11], "code": [6, 8, 9, 10, 11], "check": 6, "user": [6, 11], "guid": [6, 11], "ha": [6, 8, 9], "ocean_hgrid": 6, "addit": 6, "coordinates_north": 6, "todo": 6, "describ": 6, "what": 6, "wa": 6, "parallel": 6, "u": [6, 7, 8, 9], "equat": 6, "etc": 6, "param": 6, "f90": [6, 8, 9, 10, 11], "modif": 6, "2": [6, 8, 9, 10], "degre": 6, "resolut": 6, "orca1": 6, "also": 6, "ori": 6, "trop": 6, "A": 6, "simpl": 6, "diff": 6, "intend": 6, "compil": 6, "under": 6, "compliant": 6, "fortran": 6, "90": [6, 7, 11], "It": 6, "flag": 6, "promot": 6, "real": 6, "byte": 6, "link": 6, "librari": 6, "casper": 6, "intel": 6, "4": [6, 7, 8, 9, 10, 11], "Then": [6, 7], "cd": 6, "clean": 6, "call": 6, "tripol": 6, "ex": 6, "three": 6, "map": [7, 11], "cvmix": 7, "tidal": 7, "parameter": 7, "xesmf": 7, "xe": 7, "plot_opt": 7, "10e": 7, "src_d": 7, "altunta": 7, "mom": 7, "tx0": 7, "66v1": 7, "gen_grid_190314": 7, "tidal_tx0": 7, "energy_new": 7, "wave_dissip": 7, "collect": [7, 11], "quadmesh": [7, 11], "0x2b831feb33a0": 7, "discard": 7, "correspond": 7, "5": [7, 8, 9, 10, 11], "deg": 7, "longitud": [7, 8, 9, 10, 11], "appli": [7, 9], "nearest_s2d": 7, "regridd": 7, "slice": 7, "method": 7, "period": 7, "src_ds_fix": 7, "0x2b8328db2370": 7, "dst_grid_path": 7, "tx2_3_grid": 7, "dst_grid": 7, "renam": 7, "tlon": 7, "lon": [7, 9, 10], "tlat": 7, "lat": [7, 9, 10], "qlon": 7, "lon_b": 7, "qlat": 7, "lat_b": 7, "dst_fld": 7, "app": [7, 8, 9], "opt": [7, 8, 9], "npl": [7, 8, 9], "2022b": 7, "lib": [7, 8, 9], "python3": [7, 8, 9], "site": [7, 8, 9], "packag": [7, 8, 9], "backend": 7, "userwarn": 7, "latitud": [7, 8, 9, 10, 11], "outsid": 7, "0x2b8329192640": 7, "tmask": 7, "0x2b83291e3490": 7, "dst_d": 7, "estim": 7, "jayn": 7, "whoi": 7, "alper": 7, "m": [7, 8, 9, 10, 11], "energy_new_tx2_3_conserve_230415": 7, "NOT": 7, "allow": [7, 8, 9], "energy_new_tx2_3_conserve_230415_cdf5": 7, "load_ext": [8, 9, 10, 11], "autoreload": [8, 9, 10, 11], "os": [8, 9, 10, 11], "sy": [8, 9, 10, 11], "insert": [8, 9, 10, 11], "abspath": [8, 9, 10, 11], "date": [8, 9, 10, 11], "topo_edit_util": [8, 9, 10, 11], "inspect_topo": [8, 9, 10], "create_soc_topo_t": [8, 9, 10], "topo": [8, 9, 10, 11], "srtm15_v2": [8, 9, 10, 11], "nsub": [8, 9, 10], "sub150": [8, 9, 10, 11], "topo_src": [8, 9, 10], "srtm": [8, 9, 10, 11], "edit4": [8, 9, 10, 11], "sml1": [8, 9, 10], "0_c1": [8, 9, 10], "depth_var_in": [8, 9, 10], "d_interp": [8, 9, 10], "path_root": [8, 9, 10], "path_in": [8, 9, 10], "file_in": [8, 9, 10, 11], "dss": [8, 9, 10], "depth_var_new": [8, 9], "depth": [8, 9, 10, 11], "copi": [8, 9, 11], "deep": [8, 9], "dimens": [8, 9, 10, 11], "540": [8, 9, 10, 11], "480": [8, 9, 10, 11], "lonq": [8, 9, 10, 11], "541": [8, 9, 10, 11], "latq": [8, 9, 10, 11], "481": [8, 9, 10, 11], "285": [8, 9, 10, 11], "284": [8, 9, 10, 11], "81": [8, 9, 10, 11], "89": [8, 9, 10, 11], "287": [8, 9, 10, 11], "31": [8, 9, 10, 11], "91": [8, 9, 10, 11], "geolonb": [8, 9, 10, 11], "geolatb": [8, 9, 10, 11], "z": [8, 9, 10, 11], "float32": [8, 9, 10, 11], "ocn_frac": [8, 9, 10, 11], "d_min": [8, 9, 10, 11], "d_max": [8, 9, 10, 11], "hand_edit": [8, 9, 10, 11], "int32": [8, 9, 10, 11], "orig_mask": [8, 9, 10, 11], "statist": [8, 9, 10, 11], "creator": [8, 9, 10, 11], "frank": [8, 9, 10, 11], "bryan": [8, 9, 10, 11], "20240216": [8, 9, 10, 11], "create_model_topo": [8, 9, 10, 11], "version": [8, 9, 10, 11], "campaign": [8, 9, 10, 11], "cgd": [8, 9, 10, 11], "oc": [8, 9, 10, 11], "srtm15": [8, 9, 10, 11], "histori": [8, 9, 10, 11], "lake": [8, 9, 10], "fill": [8, 9, 10], "2024xarrai": [8, 9, 10, 11], "datasetdimens": [8, 9, 10, 11], "540lath": [8, 9, 10, 11], "480lonq": [8, 9, 10, 11], "541latq": [8, 9, 10, 11], "481coordin": [8, 9, 10, 11], "67longnam": [8, 9, 10, 11], "center": [8, 9, 10, 11], "pointsunit": [8, 9, 10, 11], "degrees_eastarrai": [8, 9, 10, 11], "666667": [8, 9, 10, 11], "333333": [8, 9, 10, 11], "86longnam": [8, 9, 10, 11], "558125": [8, 9, 10, 11], "459689": [8, 9, 10, 11], "360114": [8, 9, 10, 11], "327983": [8, 9, 10, 11], "596848": [8, 9, 10, 11], "856897": [8, 9, 10, 11], "0longnam": [8, 9, 10, 11], "corner": [8, 9, 10, 11], "91longnam": [8, 9, 10, 11], "60692": [8, 9, 10, 11], "509049": [8, 9, 10, 11], "410045": [8, 9, 10, 11], "458332": [8, 9, 10, 11], "719382": [8, 9, 10, 11], "908629": [8, 9, 10, 11], "longnam": [8, 9, 10, 11], "degrees_east": [8, 9, 10, 11], "259200": [8, 9, 10, 11], "valu": [8, 9, 10, 11], "dtype": [8, 9, 10, 11], "degrees_north": [8, 9, 10, 11], "260221": [8, 9, 10], "averag": [8, 9, 10, 11], "elev": [8, 9, 10, 11], "land": [8, 9, 10], "cellunit": [8, 9, 10, 11], "fraction": [8, 9, 10, 11], "cover": [8, 9, 10, 11], "oceanunit": [8, 9, 10, 11], "dimensionless": [8, 9, 10, 11], "tracer": [8, 9, 10, 11], "d_mean": [8, 9, 10, 11], "mean": [8, 9, 10, 11], "d_median": [8, 9, 10, 11], "median": [8, 9, 10, 11], "d2_mean": [8, 9, 10, 11], "squar": [8, 9, 10, 11], "minimum": [8, 9, 10, 11], "depthunit": [8, 9, 10], "msmooth": [8, 9, 10], "scale": [8, 9, 10], "l": [8, 9, 10], "length": [8, 9, 10, 11], "0smooth": [8, 9, 10], "function": [8, 9, 10], "lonhpandasindexpandasindex": [8, 9, 10, 11], "66666666666663": [8, 9, 10, 11], "33333333333337": [8, 9, 10, 11], "283": [8, 9, 10, 11], "33333333333326": [8, 9, 10, 11], "282": [8, 9, 10, 11], "66666666666674": [8, 9, 10, 11], "280": [8, 9, 10, 11], "66": [8, 9, 10, 11], "66666666666652": [8, 9, 10, 11], "68": [8, 9, 10, 11], "70": [8, 9, 10, 11], "lathpandasindexpandasindex": [8, 9, 10, 11], "55812493176025": [8, 9, 10, 11], "45968892243141": [8, 9, 10, 11], "36011371758804": [8, 9, 10, 11], "25938643762682": [8, 9, 10, 11], "15749406809196": [8, 9, 10, 11], "05442345864881": [8, 9, 10, 11], "95016132206449": [8, 9, 10, 11], "84469423319594": [8, 9, 10, 11], "7380086279867": [8, 9, 10, 11], "63009080247244": [8, 9, 10, 11], "87": [8, 9, 10, 11], "35130632947785": [8, 9, 10, 11], "64265130214595": [8, 9, 10, 11], "93114657007989": [8, 9, 10, 11], "21672143255606": [8, 9, 10, 11], "49930050239352": [8, 9, 10, 11], "77879697673029": [8, 9, 10, 11], "05509417487352": [8, 9, 10, 11], "32798253720107": [8, 9, 10, 11], "59684785783186": [8, 9, 10, 11], "85689740546577": [8, 9, 10, 11], "lonqpandasindexpandasindex": [8, 9, 10, 11], "latqpandasindexpandasindex": [8, 9, 10, 11], "60691972636295": [8, 9, 10, 11], "50904852679943": [8, 9, 10, 11], "4100445212757": [8, 9, 10, 11], "30989489722798": [8, 9, 10, 11], "20858670775536": [8, 9, 10, 11], "10610687059054": [8, 9, 10, 11], "00244216707667": [8, 9, 10, 11], "89757924115176": [8, 9, 10, 11], "791504598341": [8, 9, 10, 11], "68420460475764": [8, 9, 10, 11], "49767338655563": [8, 9, 10, 11], "78730978470759": [8, 9, 10, 11], "07401954419664": [8, 9, 10, 11], "3577100200491": [8, 9, 10, 11], "63826530944964": [8, 9, 10, 11], "91551547186462": [8, 9, 10, 11], "18914498663644": [8, 9, 10, 11], "4583322856567": [8, 9, 10, 11], "7193824871921": [8, 9, 10, 11], "90862896010506": [8, 9, 10, 11], "gridcreat": [8, 9, 10, 11], "20240216gener": [8, 9, 10, 11], "f90model": [8, 9, 10, 11], "tx2_3v2sourc": [8, 9, 10, 11], "ncedit": [8, 9, 10, 11], "iedit": [8, 9], "jedit": [8, 9], "zedit": [8, 9], "place": [8, 9, 10], "soc_tabl": [8, 9, 10], "lon_beg": [8, 9, 10], "236": [8, 11], "lon_end": [8, 9, 10], "231": 8, "lat_beg": [8, 9, 10], "lat_end": [8, 9, 10], "zmax": [8, 9, 10], "2000": 8, "i": [8, 9, 11], "j": [8, 9, 11], "260": 8, "900": 8, "500": [8, 9, 10], "600": [8, 9], "950": 8, "1000": [8, 9, 10], "1100": 8, "1200": 8, "float": [8, 9], "n": [8, 9, 10], "shape": [8, 9], "concaten": [8, 9], "238": [8, 9, 11], "292": 8, "293": [8, 11], "269": 8, "black": [8, 9, 10], "28": [8, 9, 10, 11], "40": [8, 9, 10, 11], "473": [8, 9, 11], "355": [8, 9, 11], "dardannel": [8, 9], "469": [8, 9, 11], "470": [8, 9, 11], "path_out": [8, 9, 10], "file_out_topo": [8, 9], "ocean_topo_": 8, "isoformat": 8, "ocean_topo_tx2_3v2_240501": 8, "nlong": [8, 9], "nlatg": [8, 9], "nedit": [8, 9], "ds_edit": 8, "astyp": [8, 9], "new": [8, 9], "without": [8, 9], "262": 8, "263": 8, "265": [8, 11], "ncxarrai": 8, "25coordin": 8, "int3278": 8, "470long_nam": 8, "deptharrai": [8, 9], "int32260": 8, "353long_nam": 8, "266": 8, "267": 8, "268": 8, "271": 8, "272": 8, "273": 8, "274": 8, "275": [8, 11], "276": [8, 11], "277": 8, "float64900": 8, "45long_nam": 8, "metersarrai": [8, 9], "6712265": 8, "45059586": 8, "editsorigin": [8, 9], "ds_out": [8, 9], "merg": 8, "481nedit": 8, "manual": 8, "By": 8, "url": 8, "http": 8, "github": 8, "com": 8, "ncar": 8, "append_topo_edit": 8, "ipynb": 8, "2024a": [8, 9], "dask": [8, 9], "config": [8, 9], "742": [8, 9], "futurewarn": [8, 9], "kei": [8, 9], "failur": [8, 9], "deprec": [8, 9], "distribut": [8, 9], "schedul": [8, 9], "depth_new": 9, "interpol": [9, 10], "todai": [9, 10, 11], "topoedit_": 9, "topoedit_tx2_3v2_02": 9, "width": 9, "around": 9, "391": 9, "408": [9, 11], "59": [9, 11], "1500": [9, 10], "800": 9, "417": 9, "416": 9, "415": [9, 11], "396": [9, 11], "395": [9, 11], "397": [9, 11], "398": [9, 11], "1300": [9, 10], "92": [9, 10, 11], "421": 9, "420": 9, "947": 9, "60": [9, 11], "309": 9, "310": [9, 11], "311": [9, 11], "320": [9, 11], "322": [9, 11], "321": [9, 11], "323": [9, 11], "325": [9, 11], "326": [9, 11], "327": 9, "328": [9, 11], "2500": 9, "1650": 9, "315": 9, "64": [9, 11], "1800": 9, "312": [9, 11], "sea": [9, 10], "242": [9, 10, 11], "350": [9, 10], "219": [9, 11], "137": [9, 10, 11], "493": 9, "496": 9, "305": [9, 11], "171": 9, "167": 9, "85": [9, 11], "387": [9, 11], "388": 9, "390": [9, 11], "494": [9, 11], "495": [9, 11], "ni": 9, "nj": 9, "assign_attr": 9, "ss": [9, 11], "303": [9, 11], "304": 9, "302": [9, 11], "522": 9, "515": 9, "512": 9, "658": 9, "251": [9, 11], "57": [9, 11], "56coordin": 9, "int32386": 9, "495long_nam": 9, "int32408": 9, "304long_nam": 9, "float64522": 9, "78long_nam": 9, "50872803": 9, "45153809": 9, "74291992": 9, "732": 9, "344": [9, 11], "06497192": 9, "372": [9, 11], "04974365": 9, "318": [9, 11], "84753418": 9, "815": 9, "05352783": 9, "460": [9, 11], "15078735": 9, "359": 9, "96356201": 9, "28036499": 9, "807": 9, "59130859": 9, "510": 9, "43533325": 9, "419": [9, 11], "18014526": 9, "56967163": 9, "671": 9, "648": 9, "05291748": 9, "608": 9, "49169922": 9, "642": 9, "29437256": 9, "1221": 9, "54528809": 9, "1199": 9, "88659668": 9, "1237": 9, "68579102": 9, "1190": 9, "57165527": 9, "1109": 9, "51916504": 9, "1073": 9, "81286621": 9, "949": 9, "30462646": 9, "868": 9, "891": 9, "5703125": 9, "1096": 9, "52990723": 9, "1268": 9, "77990723": 9, "1008": 9, "82556152": 9, "893": 9, "60858154": 9, "890": 9, "99505615": 9, "941": 9, "45379639": 9, "864": 9, "677": 9, "869": 9, "23199463": 9, "873": 9, "792": 9, "1112": 9, "830": 9, "878": 9, "1667": 9, "78100586": 9, "1669": 9, "08996582": 9, "22503662": 9, "66130066": 9, "399": 9, "96154785": 9, "9463501": 9, "18634033": 9, "45144653": 9, "78413773": 9, "int32480arrai": 9, "int32540arrai": 9, "unlimited_dim": 9, "edit": 10, "hand": [10, 11], "file_out_chan": 10, "channels_": 10, "txt": 10, "fmt_out": 10, "2f": 10, "1f": 10, "channels_tx2_3v2_02": 10, "st": 10, "lon1": 10, "lat1": 10, "lon2": 10, "lat2": 10, "xbox": 10, "ybox": 10, "linestyl": 10, "dash": 10, "line": 10, "line2d": 10, "0x147a63960210": 10, "0e3": 10, "w": 10, "u_width": 10, "write": 10, "12000": 10, "0x147a5b910650": 10, "v_width": 10, "5000": 10, "0x147a5b841090": 10, "115": 10, "0x147a5b5d3dd0": 10, "22000": 10, "0x147a5ae2b710": 10, "32000": 10, "region": 11, "shown": 11, "flip": 11, "enter": 11, "binari": 11, "algorithm": 11, "inland": 11, "organ": 11, "basin": 11, "within": 11, "notebook": 11, "map_mask": 11, "mom6_latlon2ij": 11, "mask_flood": 11, "filterwarn": 11, "ignor": 11, "file_root": 11, "topo_data": 11, "edit_no": 11, "df_old": 11, "file_out": 11, "df_new": 11, "pcolormesh": 11, "180": 11, "lonbeg": 11, "lonend": 11, "latbeg": 11, "latend": 11, "skp": 11, "441": 11, "437": 11, "385": 11, "439": 11, "383": 11, "440": 11, "381": 11, "436": 11, "379": 11, "382": 11, "435": 11, "377": 11, "380": 11, "434": 11, "376": 11, "433": 11, "po": 11, "362": 11, "447": 11, "467": 11, "341": 11, "468": 11, "456": 11, "354": 11, "370": 11, "371": 11, "331": 11, "364": 11, "366": 11, "335": 11, "363": 11, "336": 11, "316": 11, "98": 11, "99": 11, "96": 11, "175": 11, "174": 11, "352": 11, "172": 11, "153": 11, "333": 11, "497": 11, "514": 11, "516": 11, "248": 11, "234": 11, "257": 11, "161": 11, "154": 11, "195": 11, "317": 11, "193": 11, "250": 11, "255": 11, "253": 11, "296": 11, "133": 11, "134": 11, "125": 11, "237": 11, "239": 11, "47": 11, "373": 11, "244": 11, "123": 11, "113": 11, "76": 11, "224": 11, "337": 11, "339": 11, "165": 11, "166": 11, "158": 11, "set_ylim": 11, "132": 11, "126": 11, "459": 11, "206": 11, "208": 11, "461": 11, "211": 11, "464": 11, "216": 11, "222": 11, "221": 11, "104": 11, "474": 11, "229": 11, "475": 11, "228": 11, "479": 11, "225": 11, "313": 11, "477": 11, "209": 11, "105": 11, "458": 11, "463": 11, "368": 11, "455": 11, "357": 11, "378": 11, "457": 11, "374": 11, "84": 11, "451": 11, "343": 11, "426": 11, "428": 11, "365": 11, "422": 11, "452": 11, "425": 11, "394": 11, "407": 11, "406": 11, "438": 11, "462": 11, "431": 11, "476": 11, "289": 11, "478": 11, "471": 11, "446": 11, "101": 11, "448": 11, "102": 11, "143": 11, "138": 11, "ad": 11, "2023": 11, "217": 11, "207": 11, "107": 11, "110": 11, "100": 11, "155": 11, "0x150d2c284b50": 11, "0x150d2c251750": 11, "ise": 11, "jseed": 11, "seed": 11, "filled_mask": 11, "96473328390503": 11, "cpu": 11, "ms": 11, "total": 11, "241": 11, "wall": 11, "717": 11, "0x150d27623160": 11, "04": 11, "951264": 11, "909058": 11, "977755": 11, "955827": 11, "degrees_northarrai": 11, "974192": 11, "034317": 11, "182496": 11, "987185": 11, "03778": 11, "185691": 11, "038869": 11, "186724": 11, "int320": 11, "pointsarrai": 11, "0arrai": 11, "846": 11, "0x150d27b42c50": 11, "robust": 11, "0x150d2c135f00": 11}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"2": [0, 2, 11], "3": 0, "degre": 0, "global": [0, 11], "configur": 0, "contribut": 0, "reproduc": 0, "environ": 0, "introduct": [1, 5], "defin": 2, "maximum": 2, "layer": 2, "thick": 2, "h_max": 2, "z": 2, "ic": 2, "base": 2, "valu": 2, "abov": 2, "90": 2, "quantil": 2, "map": [2, 4], "function": 2, "fit": 2, "tanh": 2, "past": 2, "rms_90": 2, "below": 2, "manual": 2, "adjust": 2, "small": 2, "bottom": 2, "perform": 2, "curv": 2, "make": 2, "sure": 2, "max": 2, "first": 2, "4": 2, "5": 2, "m": 2, "save": [2, 7], "file": [2, 4, 7, 8, 9, 10], "gener": 3, "an": 3, "esmf": 3, "mesh": 3, "step": 3, "runoff": 4, "jra": 4, "55": 4, "us": 4, "c": 4, "g": [4, 11], "compset": 4, "r05": 4, "b": 4, "greenland": [4, 11], "4km": 4, "supergrid": 6, "about": 6, "usag": 6, "wave": 7, "dissip": 7, "dataset": [7, 8, 9], "origin": [7, 11], "forc": 7, "data": [7, 8, 9, 10, 11], "fix": 7, "bad": 7, "column": 7, "src": 7, "target": 7, "grid": [7, 10], "conserv": 7, "regrid": 7, "reappli": 7, "mask": [7, 11], "append": 8, "hand": [8, 9], "edit": [8, 9, 11], "topographi": [8, 9, 10], "modif": 8, "interpol": 8, "get": [8, 9, 10], "input": [8, 9, 10, 11], "set": [8, 9, 10, 11], "chang": 8, "sanhih": 8, "island": [8, 11], "between": 8, "philippin": 8, "indonesia": [8, 11], "sulu": 8, "sea": [8, 11], "bosphoru": [8, 9, 10], "dardanel": [8, 9, 10, 11], "up": [8, 9, 10], "output": [8, 9, 10, 11], "creat": [8, 9], "select": [9, 10], "sill": 9, "strait": [9, 10, 11], "denmark": 9, "faro": 9, "bank": [9, 11], "channel": [9, 10, 11], "gibralt": [9, 10, 11], "florida": [9, 11], "windward": [9, 11], "passag": [9, 11], "anegada": 9, "lombok": [9, 10, 11], "bab": [9, 10, 11], "el": [9, 10], "mandeb": [9, 10], "bere": [9, 11], "st": [9, 11], "sub": 10, "width": 10, "ocean": 11, "tx2_3v2": 11, "eastern": 11, "atlant": 11, "western": 11, "baltic": 11, "kattegat": 11, "bothnia": 11, "finland": 11, "jutland": 11, "elb": 11, "river": 11, "delta": 11, "netherland": 11, "wadden": 11, "belgium": 11, "english": 11, "namibia": 11, "mediterranean": 11, "n": 11, "adriat": 11, "aegean": 11, "cypru": 11, "bosphori": 11, "lawrenc": 11, "bai": 11, "fundi": 11, "princ": 11, "edward": 11, "long": 11, "sound": 11, "delawar": 11, "chesapeak": 11, "bahama": 11, "puerto": 11, "rico": 11, "mona": 11, "mexico": 11, "matagorda": 11, "s": 11, "padr": 11, "campech": 11, "yucatan": 11, "caribbean": 11, "beliz": 11, "hondura": 11, "boca": 11, "del": 11, "toro": 11, "venezuela": 11, "lake": 11, "maracaibo": 11, "aruba": 11, "amazon": 11, "mouth": 11, "porto": 11, "alegr": 11, "madryn": 11, "tierra": 11, "fuego": 11, "e": 11, "w": 11, "indian": 11, "red": 11, "al": 11, "mandab": 11, "northern": 11, "persian": 11, "gulf": 11, "hormuz": 11, "gang": 11, "brahmaputra": 11, "indonesian": 11, "throghflow": 11, "sw": 11, "sumba": 11, "se": 11, "timor": 11, "sulawesi": 11, "pacif": 11, "hawaii": 11, "galapago": 11, "alaska": 11, "glacier": 11, "ketchican": 11, "vancouv": 11, "juan": 11, "de": 11, "fuca": 11, "puget": 11, "columbia": 11, "baja": 11, "california": 11, "isla": 11, "cedro": 11, "costa": 11, "rica": 11, "panama": 11, "chile": 11, "corcovado": 11, "japan": 11, "kyushu": 11, "alaskan": 11, "arctic": 11, "cape": 11, "krusenstern": 11, "mackenzi": 11, "canadian": 11, "eskimo": 11, "nw": 11, "territori": 11, "paulatuk": 11, "victoria": 11, "goja": 11, "havn": 11, "axel": 11, "heiberg": 11, "ellesmer": 11, "robeson": 11, "central": 11, "baffin": 11, "labrador": 11, "nuugaatslaq": 11, "norwegian": 11, "ne": 11, "svalbard": 11, "iceland": 11, "norwai": 11, "narvik": 11, "russian": 11, "white": 11, "kara": 11, "laptev": 11, "antarctica": 11, "byrd": 11, "land": 11, "adeli": 11, "thwait": 11, "pine": 11, "west": 11, "smylei": 11, "antarct": 11, "pennisula": 11, "maud": 11, "fill": 11, "displai": 11, "final": 11}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinx": 56}})