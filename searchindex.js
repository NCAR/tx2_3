Search.setIndex({"docnames": ["README", "initial_conditions/ics", "mesh/mesh", "runoff_mapping/runoff", "salinity_restoring/sss_restoring", "supergrid/orca_gridgen", "tidal_energy_dissipation/regrid_tidal_dissipation", "tidal_energy_dissipation/tidal_dissipation", "topography/MaskEdit_tx2_3v2b", "topography/topo"], "filenames": ["README.md", "initial_conditions/ics.md", "mesh/mesh.md", "runoff_mapping/runoff.md", "salinity_restoring/sss_restoring.md", "supergrid/orca_gridgen.md", "tidal_energy_dissipation/regrid_tidal_dissipation.ipynb", "tidal_energy_dissipation/tidal_dissipation.md", "topography/MaskEdit_tx2_3v2b.ipynb", "topography/topo.md"], "titles": ["2/3 degree global configuration", "Introduction", "Generate an ESMF mesh", "Runoff mapping files", "Introduction", "Supergrid", "Wave Dissipation dataset", "Introduction", "Ocean Mask Edits for tx2_3v2", "Topography generation"], "terms": {"thi": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], "repositori": 0, "document": 0, "all": [0, 3, 8, 9], "step": [0, 3, 9], "need": [0, 2, 5, 6], "creat": [0, 2, 5, 6, 8, 9], "nomin": [0, 2, 5, 8], "us": [0, 2, 5, 9], "cesm": [0, 3], "mom6": [0, 5, 6, 8, 9], "first": [0, 6], "clone": 0, "your": [0, 1, 2, 4, 7, 9], "account": 0, "own": 0, "branch": 0, "work": [0, 3, 6, 8], "git": 0, "checkout": 0, "b": [0, 9], "nameofyourbranch": 0, "conda": [0, 6], "env": [0, 6], "f": [0, 9], "yml": 0, "activ": 0, "environment_nam": 0, "an": [1, 3, 4, 5, 7, 8, 9], "overview": [1, 4, 7], "project": [1, 4, 7], "feel": [1, 4, 7], "free": [1, 4, 7], "updat": [1, 4, 7], "accordingli": [1, 3, 4, 7], "todo": [2, 5, 9], "why": 2, "supergrid": 2, "ocean": [2, 3, 5, 9], "mask": [2, 3], "topographi": [2, 8], "netcdf": [2, 3, 5, 6, 9], "file": [2, 5, 9], "grid": [2, 3, 5, 8, 9], "model": [2, 5, 8], "python": 2, "gen_nc_grid": 2, "py": [2, 6, 9], "convert": [2, 9], "scrip": [2, 3], "convent": 2, "modul": [2, 5, 9], "load": [2, 5, 9], "ncl": 2, "gen_scrip": 2, "qsub": [2, 9], "create__mesh": 2, "If": 2, "script": [2, 9], "run": [2, 3, 8, 9], "successfulli": 2, "you": 2, "should": [2, 5, 9], "see": [2, 3, 9], "nc": [2, 3, 5, 6, 8, 9], "directori": 2, "cdf5": [2, 6], "nccopi": [2, 6], "k": [2, 6], "tx2_3_mesh_yymmdd": 2, "tx2_3_mesh_yymmdd_cdf5": 2, "In": [3, 9], "section": [3, 9], "we": [3, 5, 9], "generag": 3, "spread": 3, "liquid": 3, "frozen": 3, "water": [3, 8], "cime": 3, "provid": [3, 9], "tool": 3, "help": 3, "gener": [3, 5, 8], "pleas": [3, 5, 9], "build": 3, "runoff_map": 3, "execut": [3, 5, 9], "follow": [3, 5, 8, 9], "instruct": [3, 9], "here": 3, "befor": 3, "proceed": 3, "next": [3, 5], "exampl": [3, 9], "namelist": [3, 9], "from": [3, 6, 8, 9], "tx2_3": [3, 5], "map_jra_to_tx2_3": 3, "nml": 3, "input_nml": 3, "gridtyp": 3, "ob": [3, 9], "file_roff": 3, "glade": [3, 6, 8, 9], "p": [3, 9], "omwg_dev": 3, "jra55": 3, "domain": 3, "190212": 3, "roff": 3, "jra025m": 3, "190213": 3, "file_ocn": 3, "gmarqu": 3, "mesh": [3, 6], "tx2_3_scrip_yymmdd": 3, "file_nn": 3, "map_jra_to_tx2_3_nn": 3, "yymmdd": 3, "file_smooth": 3, "map_jra_to_tx2_3_sm_e333r100_yymmdd": 3, "file_new": 3, "map_jra_to_tx2_3_nnsm_e333r100_yymmdd": 3, "titl": 3, "nearest": 3, "neighbor": 3, "smooth": 3, "efold": 3, "1000000": 3, "0": [3, 6, 8, 9], "rmax": 3, "300000": 3, "step1": 3, "true": [3, 6, 8], "step2": 3, "step3": 3, "note": [3, 6, 9], "make": [3, 5, 6, 9], "sure": [3, 6, 9], "modifi": [3, 5, 9], "where": [3, 6, 8, 9], "variabl": [3, 8], "can": [3, 5, 9], "divid": 3, "four": 3, "categori": 3, "input": [3, 6, 9], "type": [3, 5], "rtm": 3, "720": [3, 6], "x": 3, "360": [3, 8], "ascii": 3, "xc": 3, "yc": 3, "xv": 3, "yv": 3, "area": [3, 8, 9], "name": 3, "must": [3, 5, 9], "contain": 3, "grid_area": 3, "along": 3, "typic": 3, "rdirc": 3, "OR": 3, "1": [3, 5, 6, 8, 9], "cell": [3, 8, 9], "3": [3, 5, 6, 8, 9], "below": 3, "file_ocn_coastal_mask": 3, "onli": [3, 5], "coastal": 3, "specifi": 3, "The": [3, 5, 8, 9], "file_ocn_coast_mask": 3, "standard": 3, "includ": [3, 5], "paramet": 3, "distanc": 3, "meter": [3, 8], "default": [3, 5, 8, 9], "maximum": [3, 8], "radiu": 3, "effect": 3, "set": [3, 9], "string": 3, "add": [3, 8, 9], "unset": 3, "restrict_smooth_src_to_nn_dest": 3, "option": 3, "limit": 3, "sourc": [3, 8], "point": [3, 5, 8, 9], "just": 3, "get": 3, "fals": 3, "instead": 3, "comput": 3, "multipli": 3, "two": 3, "togeth": 3, "output": 3, "field": [3, 6], "nn": 3, "smoother": 3, "combin": 3, "nnsm": 3, "jra_to_ocn": 3, "sh": 3, "sandbox": 3, "cesm2_3_beta08_carib": 3, "gen_mapping_fil": 3, "runoff_to_ocn": 3, "src": 3, "env_mach_specif": 3, "export": 3, "tmpdir": 3, "scratch": [3, 8], "map_r05_to_tx2_3": 3, "cseg": 3, "inputdata": 3, "lnd": 3, "clm2": 3, "rtmdata": 3, "05": 3, "061026": 3, "tx2_3_scrip_230415": 3, "map_r05_to_tx2_3_nn": 3, "230415": 3, "map_r05_to_tx2_3_sm_e1000r1000_230415": 3, "map_r05_to_tx2_3_nnsm_e1000r1000_230415": 3, "r05_to_ocn": 3, "gland4km": 3, "map_greenland_4km_to_tx2_3": 3, "glc": 3, "cism": 3, "griddata": 3, "scripgrid_greenland_4km_epsg3413_c170414": 3, "map_gland4km_to_tx2_3_nn": 3, "map_gland4km_to_tx2_3_sm_e1000r300_230415": 3, "map_gland4km_to_tx2_3_nnsm_e1000r300_230415": 3, "gland4km_to_ocn": 3, "orca_gridgen": 5, "which": 5, "reli": 5, "nemo": 5, "framework": 5, "tripolar": [5, 8], "For": [5, 8, 9], "complet": 5, "descript": [5, 6, 8], "code": [5, 8], "check": 5, "user": 5, "guid": [5, 8], "origin": 5, "ha": 5, "been": 5, "ocean_hgrid": 5, "addit": 5, "coordin": [5, 8], "coordinates_north": 5, "describ": 5, "what": 5, "els": 5, "wa": 5, "parallel": 5, "u": [5, 6], "equat": 5, "etc": 5, "param": 5, "f90": [5, 8, 9], "modif": [5, 9], "2": [5, 9], "degre": 5, "resolut": [5, 9], "configur": 5, "orca1": 5, "also": 5, "comparison": 5, "ori": 5, "trop": 5, "A": 5, "simpl": 5, "diff": 5, "chang": [5, 9], "intend": 5, "compil": 5, "under": 5, "ani": [5, 8], "compliant": 5, "fortran": 5, "90": [5, 6, 8], "It": 5, "flag": 5, "promot": 5, "real": 5, "8": [5, 6, 8, 9], "byte": 5, "link": 5, "librari": 5, "To": [5, 9], "casper": [5, 9], "intel": 5, "19": [5, 8, 9], "4": [5, 6, 8, 9], "Then": [5, 6], "cd": 5, "clean": 5, "call": 5, "tripol": 5, "ex": 5, "three": 5, "map": [6, 8], "cvmix": 6, "tidal": 6, "parameter": 6, "import": [6, 8], "xarrai": [6, 8], "xr": [6, 8], "xesmf": 6, "xe": 6, "numpi": [6, 8], "np": [6, 8], "matplotlib": [6, 8], "color": 6, "datetim": [6, 8], "plot_opt": 6, "size": 6, "norm": 6, "lognorm": 6, "vmin": [6, 8], "10e": 6, "7": [6, 8], "vmax": [6, 8], "src_d": 6, "open_dataset": [6, 8], "altunta": 6, "mom": 6, "tx0": 6, "66v1": 6, "gen_grid_190314": 6, "tidal_tx0": 6, "energy_new": 6, "wave_dissip": 6, "plot": [6, 8], "collect": [6, 8], "quadmesh": [6, 8], "0x2ac37932e850": 6, "discard": 6, "correspond": 6, "5": [6, 8, 9], "deg": 6, "longitud": [6, 8], "appli": [6, 9], "nearest_s2d": 6, "regridd": 6, "isel": 6, "slice": 6, "method": 6, "period": 6, "src_ds_fix": 6, "0x2ac38223c1c0": 6, "dst_grid_path": 6, "tx2_3_grid": 6, "dst_grid": 6, "renam": 6, "tlon": 6, "lon": 6, "tlat": 6, "lat": 6, "qlon": 6, "lon_b": 6, "qlat": 6, "lat_b": 6, "dst_fld": 6, "app": 6, "opt": 6, "npl": 6, "2022b": 6, "lib": 6, "python3": 6, "site": 6, "packag": 6, "backend": 6, "53": [6, 8], "userwarn": 6, "latitud": [6, 8], "outsid": 6, "warn": [6, 8], "0x2ac38258bbb0": 6, "tmask": 6, "0x2ac3823d6b80": 6, "dst_d": 6, "coord": 6, "attr": [6, 8], "estim": 6, "s": 6, "r": [6, 9], "jayn": 6, "whoi": 6, "author": 6, "alper": 6, "ucar": 6, "edu": 6, "now": 6, "strftime": [6, 8], "y": [6, 8], "d": [6, 8, 9], "m": [6, 8], "h": [6, 9], "to_netcdf": [6, 8], "energy_new_tx2_3_conserve_011023": 6, "format": [6, 8], "NOT": 6, "allow": [6, 8], "energy_new_tx2_3_conserve_011023_cdf5": 6, "each": 8, "region": 8, "fraction": [8, 9], "shown": 8, "i": 8, "j": [8, 9], "label": 8, "flip": 8, "ar": [8, 9], "enter": 8, "hand": [8, 9], "binari": 8, "algorithm": 8, "inland": 8, "organ": 8, "basin": 8, "within": 8, "notebook": [8, 9], "load_ext": 8, "autoreload": 8, "inlin": 8, "os": 8, "date": 8, "pyplot": 8, "plt": 8, "map_mask": [8, 9], "mom6_latlon2ij": 8, "filterwarn": 8, "ignor": 8, "file_root": 8, "topo": 8, "sub150": 8, "topo_data": 8, "srtm": [8, 9], "edit_no": 8, "file_in": 8, "path": 8, "bryan": 8, "print": 8, "df_old": 8, "file_out": 8, "edit4": 8, "df_new": 8, "copi": 8, "def": 8, "ice9": 8, "xcyclic": 8, "iter": 8, "stack": 8, "base": 8, "implement": 8, "ic": 8, "9": 8, "flood": 8, "start": 8, "treat": 8, "posit": 8, "valu": 8, "passabl": 8, "zero": 8, "neg": 8, "block": 8, "cyclic": 8, "behavior": 8, "last": 8, "index": 8, "fold": 8, "across": 8, "top": 8, "most": 8, "edg": 8, "return": 8, "arrai": 8, "wetmask": 8, "nj": 8, "ni": 8, "shape": 8, "while": 8, "pop": 8, "continu": 8, "elif": 8, "tri": 8, "polar": 8, "fig": 8, "ax": 8, "subplot": 8, "figsiz": 8, "12": 8, "pcolormesh": 8, "geolonb": 8, "geolatb": 8, "set_xlim": 8, "lonh": 8, "180": 8, "lonbeg": 8, "lonend": 8, "13": 8, "latbeg": 8, "latend": 8, "58": 8, "skp": 8, "ocn_frac": 8, "387": 8, "441": 8, "18": 8, "23": 8, "57": 8, "62": 8, "22": 8, "30": 8, "61": 8, "11": 8, "55": 8, "390": 8, "437": 8, "52": 8, "56": 8, "385": 8, "439": 8, "383": 8, "440": 8, "50": 8, "54": 8, "381": 8, "436": 8, "379": 8, "382": 8, "435": 8, "377": 8, "380": 8, "434": 8, "48": 8, "376": 8, "433": 8, "34": 8, "38": 8, "skip": 8, "17": 8, "26": 8, "14": 8, "43": 8, "46": 8, "po": 8, "362": 8, "447": 8, "41": 8, "342": 8, "467": 8, "341": 8, "468": 8, "353": 8, "456": 8, "25": 8, "39": 8, "42": 8, "355": 8, "473": 8, "354": 8, "470": 8, "70": 8, "51": 8, "370": 8, "328": 8, "371": 8, "329": 8, "331": 8, "66": 8, "60": 8, "364": 8, "366": 8, "335": 8, "363": 8, "336": 8, "74": 8, "69": 8, "40": 8, "325": 8, "321": 8, "323": 8, "77": 8, "36": 8, "344": 8, "348": 8, "316": 8, "82": 8, "72": 8, "20": 8, "28": 8, "75": 8, "65": 8, "312": 8, "320": 8, "322": 8, "98": 8, "95": 8, "29": 8, "326": 8, "284": 8, "99": 8, "96": 8, "27": 8, "92": 8, "88": 8, "318": 8, "303": 8, "89": 8, "85": 8, "15": [8, 9], "83": 8, "16": 8, "81": 8, "10": [8, 9], "73": 8, "68": 8, "49": 8, "175": 8, "174": 8, "352": 8, "172": 8, "351": 8, "63": 8, "44": 8, "153": 8, "333": 8, "67": 8, "302": 8, "494": 8, "497": 8, "300": 8, "305": 8, "495": 8, "32": 8, "324": 8, "514": 8, "516": 8, "270": 8, "265": 8, "248": 8, "242": 8, "219": 8, "64": 8, "223": 8, "243": 8, "238": 8, "6": [8, 9], "233": 8, "234": 8, "251": 8, "257": 8, "236": 8, "161": 8, "154": 8, "314": 8, "196": 8, "195": 8, "317": 8, "193": 8, "191": 8, "93": 8, "250": 8, "255": 8, "293": 8, "294": 8, "253": 8, "296": 8, "137": 8, "133": 8, "59": 8, "134": 8, "130": 8, "129": 8, "125": 8, "237": 8, "239": 8, "121": 8, "47": 8, "373": 8, "244": 8, "123": 8, "45": 8, "116": 8, "113": 8, "86": 8, "78": 8, "76": 8, "71": 8, "224": 8, "337": 8, "339": 8, "170": 8, "165": 8, "166": 8, "158": 8, "set_ylim": 8, "87": 8, "132": 8, "126": 8, "459": 8, "206": 8, "460": 8, "208": 8, "461": 8, "211": 8, "464": 8, "216": 8, "222": 8, "221": 8, "111": 8, "104": 8, "474": 8, "229": 8, "475": 8, "228": 8, "479": 8, "225": 8, "313": 8, "310": 8, "477": 8, "209": 8, "105": 8, "458": 8, "463": 8, "91": 8, "469": 8, "368": 8, "79": 8, "80": 8, "455": 8, "357": 8, "378": 8, "457": 8, "372": 8, "374": 8, "84": 8, "451": 8, "413": 8, "415": 8, "343": 8, "409": 8, "426": 8, "428": 8, "365": 8, "422": 8, "396": 8, "452": 8, "397": 8, "425": 8, "394": 8, "419": 8, "407": 8, "408": 8, "398": 8, "406": 8, "395": 8, "ncol": 8, "287": 8, "438": 8, "462": 8, "431": 8, "476": 8, "289": 8, "275": 8, "478": 8, "471": 8, "245": 8, "446": 8, "101": 8, "448": 8, "102": 8, "143": 8, "138": 8, "ad": 8, "2023": 8, "217": 8, "207": 8, "107": 8, "117": 8, "110": 8, "100": 8, "276": 8, "311": 8, "334": 8, "345": 8, "103": 8, "338": 8, "sum": 8, "hand_edit": 8, "155": 8, "0x2b32f206f2e0": 8, "ise": 8, "jseed": 8, "seed": 8, "geolon": 8, "geolat": 8, "filled_mask": 8, "283": 8, "33333333333326": 8, "9998646381081644": 8, "orig_mask": 8, "todai": 8, "ss": 8, "histori": 8, "04": 8, "lt": 8, "dataset": [8, 9], "gt": 8, "dimens": 8, "540": 8, "lath": 8, "480": 8, "lonq": 8, "541": 8, "latq": 8, "481": 8, "float64": 8, "286": 8, "285": 8, "33": 8, "31": 8, "z": 8, "float32": 8, "d_median": 8, "d2_mean": 8, "d_min": 8, "d_max": 8, "int32": 8, "attribut": 8, "statist": [8, 9], "creator": 8, "frank": 8, "20230405": 8, "create_model_topo": [8, 9], "version": 8, "observ": 8, "2023xarrai": 8, "datasetdimens": 8, "540lath": 8, "480lonq": 8, "541latq": 8, "481coordin": 8, "67longnam": 8, "center": 8, "pointsunit": 8, "degrees_eastarrai": 8, "666667": 8, "333333": 8, "86longnam": 8, "558125": 8, "459689": 8, "360114": 8, "327983": 8, "596848": 8, "856897": 8, "0longnam": 8, "corner": 8, "91longnam": 8, "60692": 8, "509049": 8, "410045": 8, "458332": 8, "719382": 8, "908629": 8, "longnam": 8, "degrees_east": 8, "259200": 8, "dtype": 8, "degrees_north": 8, "951264": 8, "909058": 8, "977755": 8, "955827": 8, "degrees_northarrai": 8, "974192": 8, "034317": 8, "182496": 8, "987185": 8, "03778": 8, "185691": 8, "038869": 8, "186724": 8, "averag": 8, "elev": 8, "cellunit": 8, "cover": [8, 9], "oceanunit": 8, "dimensionless": 8, "int320": 8, "tracer": 8, "pointsarrai": 8, "d_mean": 8, "mean": [8, 9], "depth": 8, "median": [8, 9], "squar": [8, 9], "minimum": [8, 9], "0arrai": 8, "gridcreat": 8, "bryancreat": 8, "20230405gener": 8, "f90model": 8, "tx2_3v2sourc": 8, "srtm15_v2": [8, 9], "ncedit": 8, "846": 8, "0x2b32f21666d0": 8, "robust": 8, "0x2b32f22178e0": 8, "walk": 9, "through": 9, "requir": 9, "accur": 9, "represent": 9, "These": 9, "summar": 9, "visuali": 9, "inspect": 9, "manual": 9, "via": 9, "jupyt": 9, "high": 9, "global": 9, "interp_smooth": 9, "coupl": 9, "util": 9, "constant": 9, "kind": 9, "mom6_grid": 9, "ncdf_wrapper": 9, "gnumakefil": 9, "expect": 9, "gnu": 9, "On": 9, "pb": 9, "run_create_topo": 9, "tx1_4v2": 9, "run_intep_smooth_tx1_4v2": 9, "program": 9, "edit": 9, "maskedit_tx2_3v2b": 9, "ipynb": 9, "bathymetri": 9, "arc": 9, "sec": 9, "deriv": 9, "shuttl": 9, "radar": 9, "mission": 9, "found": 9, "campaign": 9, "cgd": 9, "oc": 9, "more": 9, "inform": 9, "its": 9, "associ": 9, "public": 9, "tozer": 9, "sandwel": 9, "t": 9, "smith": 9, "w": 9, "olson": 9, "c": 9, "beal": 9, "wessel": 9, "2019": 9, "srtm15": 9, "earth": 9, "space": 9, "scienc": 9, "1847": 9, "http": 9, "doi": 9, "org": 9, "1029": 9, "2019ea000658": 9, "ncarenv": 9, "ncarcompil": 9, "mpt": 9, "gmake": 9, "inter_smooth": 9, "maximun": 9, "tx2_3v2": 9, "open": 9, "correct": 9, "sub": 9, "nsub": 9, "defin": 9, "consist": 9, "possibl": 9, "do": 9, "so": 9, "done": 9, "one": 9, "want": 9, "around": 9, "antarctica": 9, "navig": 9, "necessari": 9, "about": 9, "made": 9, "previou": 9, "run_intep_smooth_tx2_3v2": 9, "gerat": 9, "topodata": 9, "sfnc": 9, "final": 9, "how": 9}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"2": [0, 8], "3": 0, "degre": 0, "global": [0, 8], "configur": 0, "contribut": 0, "reproduc": 0, "environ": 0, "introduct": [1, 4, 7], "gener": [2, 9], "an": 2, "esmf": 2, "mesh": 2, "about": [2, 5], "step": 2, "runoff": 3, "map": 3, "file": [3, 6], "jra": 3, "55": 3, "us": 3, "c": 3, "g": [3, 8], "compset": 3, "r05": 3, "b": 3, "greenland": [3, 8], "4km": 3, "supergrid": 5, "usag": 5, "wave": 6, "dissip": 6, "dataset": 6, "origin": [6, 8], "forc": 6, "data": [6, 8], "fix": 6, "bad": 6, "column": 6, "src": 6, "target": 6, "grid": 6, "conserv": 6, "regrid": 6, "reappli": 6, "mask": [6, 8, 9], "save": 6, "ocean": 8, "edit": 8, "tx2_3v2": 8, "input": 8, "set": 8, "eastern": 8, "atlant": 8, "western": 8, "baltic": 8, "kattegat": 8, "bothnia": 8, "finland": 8, "jutland": 8, "elb": 8, "river": 8, "delta": 8, "netherland": 8, "wadden": 8, "sea": [8, 9], "belgium": 8, "english": 8, "channel": 8, "gibralt": 8, "namibia": 8, "mediterranean": 8, "n": 8, "adriat": 8, "aegean": 8, "cypru": 8, "bosphori": 8, "dardanel": 8, "st": 8, "lawrenc": 8, "bai": 8, "fundi": 8, "princ": 8, "edward": 8, "island": 8, "long": 8, "sound": 8, "delawar": 8, "chesapeak": 8, "florida": 8, "strait": 8, "bahama": 8, "puerto": 8, "rico": 8, "windward": 8, "mona": 8, "passag": 8, "mexico": 8, "matagorda": 8, "s": 8, "padr": 8, "campech": 8, "yucatan": 8, "caribbean": 8, "beliz": 8, "hondura": 8, "boca": 8, "del": 8, "toro": 8, "venezuela": 8, "lake": 8, "maracaibo": 8, "aruba": 8, "amazon": 8, "mouth": 8, "porto": 8, "alegr": 8, "madryn": 8, "tierra": 8, "fuego": 8, "e": 8, "w": 8, "indian": 8, "red": 8, "bab": 8, "al": 8, "mandab": 8, "northern": 8, "persian": 8, "gulf": 8, "hormuz": 8, "gang": 8, "brahmaputra": 8, "indonesian": 8, "throghflow": 8, "indonesia": 8, "sw": 8, "lombok": 8, "sumba": 8, "se": 8, "timor": 8, "sulawesi": 8, "pacif": 8, "hawaii": 8, "galapago": 8, "alaska": 8, "glacier": 8, "ketchican": 8, "vancouv": 8, "juan": 8, "de": 8, "fuca": 8, "puget": 8, "columbia": 8, "baja": 8, "california": 8, "isla": 8, "cedro": 8, "costa": 8, "rica": 8, "panama": 8, "chile": 8, "corcovado": 8, "japan": 8, "kyushu": 8, "alaskan": 8, "arctic": 8, "bere": 8, "cape": 8, "krusenstern": 8, "mackenzi": 8, "canadian": 8, "eskimo": 8, "nw": 8, "territori": 8, "paulatuk": 8, "victoria": 8, "bank": 8, "goja": 8, "havn": 8, "axel": 8, "heiberg": 8, "ellesmer": 8, "robeson": 8, "central": 8, "baffin": 8, "labrador": 8, "nuugaatslaq": 8, "norwegian": 8, "ne": 8, "svalbard": 8, "iceland": 8, "norwai": 8, "narvik": 8, "russian": 8, "white": 8, "kara": 8, "laptev": 8, "antarctica": 8, "byrd": 8, "land": [8, 9], "adeli": 8, "thwait": 8, "pine": 8, "west": 8, "smylei": 8, "antarct": 8, "pennisula": 8, "maud": 8, "fill": 8, "output": 8, "displai": 8, "final": 8, "topographi": 9, "compil": 9, "fortran": 9, "code": 9, "topo": 9, "stat": 9, "check": 9, "smooth": 9}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinx": 56}})