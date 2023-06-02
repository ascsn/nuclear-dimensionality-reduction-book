Search.setIndex({"docnames": ["HO_Time/HO_Time_Dependent", "contributors", "dft/dft", "gp/gp", "greedy/greedy", "ho/1dHO", "ho/ho", "ho/ho-rbm", "ho/ho-ritz", "introduction/emulators", "introduction/introduction", "introduction/redundancy", "landing", "markdown-notebooks", "notebooks", "rbm/lagrange", "rbm/pod", "rbm/projecting", "rbm/rbm", "rbm/training", "scattering/eim", "scattering/scattering", "transforms/transforms"], "filenames": ["HO_Time/HO_Time_Dependent.ipynb", "contributors.md", "dft/dft.md", "gp/gp.md", "greedy/greedy.md", "ho/1dHO.ipynb", "ho/ho.md", "ho/ho-rbm.md", "ho/ho-ritz.md", "introduction/emulators.md", "introduction/introduction.ipynb", "introduction/redundancy.md", "landing.md", "markdown-notebooks.md", "notebooks.ipynb", "rbm/lagrange.md", "rbm/pod.md", "rbm/projecting.md", "rbm/rbm.md", "rbm/training.md", "scattering/eim.ipynb", "scattering/scattering.ipynb", "transforms/transforms.md"], "titles": ["Application 4: Time Dependent Systems (evolution in the reduced space)", "Contributors", "Application 4: Nuclear DFT", "Application 2: The Gross-Pitaevskii Equation", "The Greedy Algorithm", "Application 1: The Quantum Harmonic Oscillator", "Application 1: The Quantum Harmonic Oscillator", "The RBM Approach", "The Ritz Variational Approach", "Why emulators?", "Introduction", "Redundancy in Simulations", "Introduction to Reduced-Basis Methods in Nuclear Physics", "Notebooks with MyST Markdown", "Content with notebooks", "Lagrange Reduced Basis", "PCA and POD", "Finding the Coefficients: Projecting", "The Reduced Basis Method", "Choosing the Basis: Training", "Application 3: The Empirical Interpolation Method", "Application 2: Two body single channel nuclear scattering", "Modifying the Basis: Stretching and Translating"], "terms": {"contribut": [0, 1, 5, 9, 10, 20, 21], "pablo": [0, 1, 9, 10, 20, 21], "giuliani": [0, 1, 9, 10, 20, 21], "kyle": [0, 1, 9, 10, 20], "godbei": [0, 1, 9, 10], "edgard": [0, 1, 10], "bonilla": [0, 1, 10], "In": [0, 9, 10, 20, 21], "thi": [0, 1, 5, 9, 10, 12, 13, 14, 20, 21], "section": [0, 9, 20, 21], "we": [0, 5, 9, 10, 20, 21], "ar": [0, 5, 9, 10, 13, 20, 21], "go": [0, 9, 20, 21], "take": [0, 5, 9, 20, 21], "look": [0, 9, 20, 21], "quantum": [0, 10, 12, 21], "us": [0, 5, 9, 10, 13, 20, 21], "here": [0, 5, 9, 12, 14, 20, 21], "an": [0, 5, 9, 10, 20, 21], "exampl": [0, 9, 10, 12, 14, 20, 21], "practic": 0, "more": [0, 10, 13, 14, 20], "complex": [0, 9], "non": [0, 20, 21], "focu": [0, 20], "harmon": [0, 12, 21], "oscil": [0, 12, 21], "our": [0, 5, 9, 10, 20, 21], "first": [0, 5, 10, 20, 21], "explor": [0, 9, 20, 21], "try": [0, 10, 20], "self": 0, "interact": [0, 9, 14, 20, 21], "particl": [0, 9], "next": [0, 9, 20, 21], "chapter": [0, 9, 20, 21], "about": [0, 9, 10, 13, 14], "linear": [0, 20, 21], "hamiltonian": [0, 5, 20, 21], "begin": [0, 5, 14, 20, 21], "equat": [0, 5, 10, 20, 21], "h": [0, 5, 21], "frac": [0, 5, 20, 21], "partial": [0, 10], "2": [0, 5, 12, 13, 20], "x": [0, 5, 9, 10, 20], "end": [0, 5, 9, 14, 20, 21], "let": [0, 9, 13, 20, 21], "drop": [0, 20, 21], "numer": [0, 5, 21], "constant": [0, 10], "simplifi": [0, 20], "notat": [0, 10, 20], "hbar": [0, 5], "1": [0, 9, 12, 14, 20, 21], "The": [0, 9, 10, 12, 13, 21], "mechan": [0, 10], "wave": [0, 20, 21], "function": [0, 5, 9, 10, 20, 21], "describ": [0, 10], "through": [0, 9, 21], "schroding": [0, 5], "left": [0, 5, 20, 21], "phi": [0, 5, 20, 21], "t": [0, 5, 10, 14, 20, 21], "right": [0, 5, 10, 12, 20, 21], "rangl": [0, 5, 20, 21], "i": [0, 5, 9, 20, 21], "over": [0, 5, 10, 20], "which": [0, 9, 10, 13, 20, 21], "result": [0, 5, 12, 20, 21], "oper": [0, 5, 20, 21], "from": [0, 5, 9, 10, 14, 20, 21], "initi": [0, 5, 20, 21], "state": [0, 5, 14], "e": [0, 5, 9, 10, 20, 21], "iht": 0, "0": [0, 5, 14, 20, 21], "To": [0, 5, 10, 20, 21], "solv": [0, 20, 21], "problem": [0, 5, 12, 21], "usual": [0, 5, 9, 10, 20], "expand": [0, 20], "its": [0, 9], "taylor": 0, "approxim": [0, 5, 9, 20, 21], "finit": [0, 5], "step": [0, 10], "delta": 0, "taken": [0, 5], "approx": [0, 20, 21], "big": [0, 5, 20], "ih": 0, "still": [0, 20], "studi": [0, 9, 10, 20], "dynam": [0, 10], "build": [0, 10, 12, 21], "transform": [0, 10], "hat": [0, 5, 9, 20, 21], "sum_k": [0, 20, 21], "n": [0, 5, 10, 14, 20, 21], "a_k": [0, 5, 20, 21], "phi_k": [0, 20, 21], "where": [0, 5, 9, 20, 21], "now": [0, 5, 9, 10, 12, 20, 21], "coeffici": [0, 20, 21], "explicit": 0, "well": [0, 5, 9, 14, 20], "ani": [0, 5, 9, 10, 13, 20, 21], "paramet": [0, 5, 9, 20, 21], "version": 0, "when": [0, 13, 20, 21], "substitut": 0, "do": [0, 5, 14, 20, 21], "galerkin": [0, 5, 10, 20, 21], "project": [0, 5, 10, 20, 21], "equaiton": [0, 21], "previou": [0, 20, 21], "have": [0, 9, 10, 13, 20, 21], "f_": [0, 5, 20, 21], "alpha": [0, 5, 9, 20, 21], "can": [0, 5, 9, 10, 12, 13, 14, 20, 21], "substitud": 0, "langl": [0, 5, 20, 21], "phi_j": [0, 20], "text": [0, 13, 20], "j": [0, 5, 10, 20, 21], "choos": [0, 5, 20], "pricinp": 0, "compon": [0, 5, 9, 20, 21], "sever": [0, 9, 10, 20], "snapshot": [0, 20], "mean": [0, 14], "delta_j": 0, "k": [0, 5, 10, 20, 21], "lead": [0, 5], "d": [0, 1, 5, 10, 13, 20, 21], "dt": 0, "a_j": 0, "boldsymbol": 0, "a_1": [0, 20, 21], "a_n": [0, 20, 21], "comput": [0, 9, 10, 20], "h_": 0, "small": [0, 20], "number": [0, 5, 20, 21], "solut": [0, 5, 20, 21], "much": [0, 9], "effici": 0, "thank": [0, 1, 10], "both": 0, "fact": 0, "ll": 0, "being": [0, 9], "dimens": [0, 20], "tradit": 0, "appreci": [0, 9], "bigger": 0, "recov": 0, "stabl": 0, "import": [0, 5, 9, 10, 14, 20, 21], "numpi": [0, 5, 14, 20, 21], "np": [0, 5, 14, 20, 21], "matplotlib": [0, 5, 14, 20, 21], "pyplot": [0, 5, 14, 20, 21], "plt": [0, 5, 14, 20, 21], "anim": 0, "funcanim": 0, "scipi": [0, 5, 20, 21], "integr": [0, 5, 20, 21], "solve_ivp": [0, 20, 21], "linalg": [0, 5, 20, 21], "expm": 0, "lstsq": 0, "interpol": [0, 12, 21], "cubicsplin": 0, "math": [0, 5, 14], "timeit": 0, "bege": 0, "creat": [0, 14, 20, 21], "deriv": [0, 20, 21], "element": [0, 5, 9, 20], "def": [0, 5, 20, 21], "generate_second_derivative_matrix": 0, "dx": [0, 5], "gener": [0, 9, 12, 20, 21], "matrix": [0, 5, 20, 21], "second": [0, 21], "five": 0, "point": [0, 5, 9, 10, 20, 21], "stencil": 0, "main_diag": 0, "ones": 0, "5": [0, 5, 9, 14, 20, 21], "off_diag": [0, 5], "3": [0, 5, 9, 12, 21], "off_diag2": 0, "12": [0, 5, 20, 21], "d2": [0, 21], "diag": [0, 5, 21], "modifi": 0, "dirichlet": 0, "boundari": [0, 5, 21], "condit": [0, 5, 20, 21], "15": [0, 5, 20], "16": [0, 21], "return": [0, 5, 20, 21], "delai": 0, "movi": 0, "sinc": [0, 5, 10, 20, 21], "produc": [0, 9, 20], "lot": [0, 14], "them": [0, 9, 10, 20, 21], "frame_per_second": 0, "60": 0, "x_max": 0, "10": [0, 5, 14, 20, 21], "maximum": 0, "coordin": [0, 20], "grid": [0, 5, 20], "n_grid": 0, "150": 0, "set": [0, 5, 9, 20], "up": [0, 5, 10], "linspac": [0, 14, 20, 21], "omega_val0": 0, "you": [0, 1, 9, 10, 12, 13, 14, 21], "chanc": [0, 9], "see": [0, 5, 10, 12, 13, 14, 20, 21], "what": [0, 10, 20, 21], "happen": [0, 9], "initial_st": 0, "exp": [0, 5, 20, 21], "sqrt": [0, 5, 20, 21], "pi": [0, 5], "harmonic_oscillator_time_evolut": 0, "total_tim": 0, "time_step": 0, "total_ord": 0, "time_step_s": 0, "potenti": [0, 5, 20, 21], "energi": [0, 5, 9, 10, 20, 21], "term": [0, 5, 20, 21], "v": [0, 5, 20, 21], "state_evolut": 0, "rho_evolut": 0, "current_st": 0, "dot": [0, 5], "conj": [0, 20], "append": 0, "real": [0, 20], "_": [0, 5, 20], "rang": [0, 5, 10, 14, 20, 21], "next_stat": 0, "copi": [0, 5, 20, 21], "current_correct": 0, "1j": 0, "arrai": [0, 5, 14, 20, 21], "hf_results0": 0, "10000": 0, "cpu": 0, "user": 0, "87": 0, "s": [0, 5, 9, 10, 13, 14, 20, 21], "sy": 0, "total": [0, 9], "6": [0, 5, 20, 21], "37": 0, "wall": 0, "seconds_tot": 0, "how": [0, 5, 9, 13, 20, 21], "long": [0, 9], "want": [0, 5, 9, 10, 14, 20], "last": [0, 10, 21], "rho_hf0": 0, "ceil": 0, "len": [0, 5, 21], "wf_hf0": 0, "movie_mak": 0, "wf_list": 0, "rholist": 0, "omega_valu": 0, "name": [0, 10, 20], "fig": [0, 5, 14, 20, 21], "figur": [0, 9, 20, 21], "dpi": [0, 20, 21], "200": [0, 20, 21], "ax": [0, 5, 14, 20, 21], "xlim": 0, "ylim": 0, "plot": [0, 5, 14, 20, 21], "color": [0, 14], "line": [0, 13, 14], "lw": [0, 14], "b": [0, 10, 20, 21], "line2": 0, "linestyl": [0, 20, 21], "dash": 0, "r": [0, 5, 10, 20, 21], "line3": 0, "orang": [0, 9], "init": [0, 13], "set_data": 0, "imag": [0, 14], "init_func": 0, "frame": 0, "interv": 0, "frame_delai": 0, "1000": 0, "blit": 0, "true": [0, 20, 21], "ffwriter": 0, "ffmpegwrit": 0, "fp": 0, "save": 0, "mp4": 0, "writer": 0, "If": [0, 1, 9, 10, 13, 20], "run": [0, 9, 13], "notebook": 0, "local": 0, "follow": [0, 13, 20, 21], "command": [0, 13], "imaginari": 0, "part": [0, 9, 10, 20, 21], "togeth": [0, 20], "densiti": [0, 10], "all": [0, 9, 13, 20, 21], "across": [0, 21], "unit": 0, "hf_solver": 0, "filenotfounderror": 0, "traceback": 0, "most": [0, 9, 20], "recent": [0, 9], "call": [0, 10, 20], "eval": [0, 5], "modul": 0, "tmp": [0, 20], "ipykernel_1774": 0, "2653633946": 0, "py": 0, "29": 0, "30": [0, 21], "31": 0, "opt": 0, "hostedtoolcach": 0, "python": [0, 10, 20], "7": 0, "x64": 0, "lib": 0, "python3": 0, "site": 0, "packag": 0, "filenam": 0, "codec": 0, "bitrat": 0, "extra_arg": 0, "metadata": 0, "extra_anim": 0, "savefig_kwarg": 0, "progress_callback": 0, "1070": 0, "widget": 0, "likewis": 0, "done": [0, 20], "savefig": 0, "1071": 0, "mpl": 0, "rc_context": 0, "bbox": 0, "none": [0, 20], "1072": 0, "_fig": 0, "1073": 0, "cbook": 0, "_setattr_cm": 0, "canva": 0, "1074": 0, "_is_sav": 0, "manag": 0, "contextlib": 0, "__enter__": 0, "110": 0, "del": 0, "arg": [0, 20, 21], "kwd": 0, "func": 0, "111": 0, "112": 0, "gen": 0, "113": 0, "except": 0, "stopiter": 0, "114": 0, "rais": 0, "runtimeerror": 0, "didn": 0, "yield": 0, "outfil": 0, "kwarg": 0, "230": 0, "231": 0, "particular": 0, "sequenc": 0, "contextmanag": 0, "232": 0, "setup": [0, 5, 9], "233": 0, "234": 0, "319": 0, "so": [0, 9, 10, 13, 20, 21], "grab_fram": 0, "write": [0, 13], "data": [0, 9, 10, 14], "pipe": 0, "320": 0, "elimin": [0, 20], "need": [0, 5, 9, 13, 20], "temp": 0, "file": [0, 13], "321": 0, "_run": 0, "322": 0, "323": 0, "331": 0, "_proc": 0, "subprocess": 0, "popen": 0, "332": 0, "stdin": 0, "stdout": 0, "stderr": 0, "333": 0, "creationflag": 0, "subprocess_creation_flag": 0, "334": 0, "335": 0, "finish": [0, 9], "__init__": 0, "bufsiz": 0, "execut": [0, 13], "preexec_fn": 0, "close_fd": 0, "shell": 0, "cwd": 0, "env": 0, "universal_newlin": 0, "startupinfo": 0, "restore_sign": 0, "start_new_sess": 0, "pass_fd": 0, "encod": 0, "error": [0, 5, 21], "798": 0, "c2pread": 0, "c2pwrite": 0, "799": 0, "errread": 0, "errwrit": 0, "800": 0, "801": 0, "802": 0, "cleanup": 0, "child": 0, "fail": [0, 20], "start": [0, 5, 9, 13], "_execute_child": 0, "p2cread": 0, "p2cwrite": 0, "1549": 0, "errno_num": 0, "errno": 0, "enoent": 0, "1550": 0, "err_msg": 0, "repr": 0, "err_filenam": 0, "1551": 0, "child_exception_typ": 0, "1552": 0, "1553": 0, "No": 0, "directori": 0, "ffmpeg": 0, "label": [0, 5, 20, 21], "g": [0, 10, 20], "linewidth": 0, "intermedi": 0, "ij": [0, 20], "55": 0, "grei": 0, "phi_0": [0, 20, 21], "final": [0, 9, 20], "t_f": 0, "legend": [0, 5, 14, 20, 21], "fontsiz": 0, "frameon": 0, "edgecolor": 0, "black": [0, 9], "xlabel": [0, 20], "20": [0, 20, 21], "show": [0, 5, 9, 13, 20, 21], "after": [0, 9], "some": [0, 1, 9, 14, 20, 21], "highlight": [0, 20], "later": [0, 21], "accuraci": [0, 9, 20], "emul": [0, 10, 20, 21], "live": [0, 1], "low": [0, 9, 20], "dimension": [0, 9, 10, 12, 20], "manifold": [0, 5], "differ": [0, 5, 9, 20, 21], "find": [0, 5, 9, 10, 12, 20], "find_best_linear_combin": 0, "y": 0, "basis_funct": 0, "column_stack": 0, "basis_func": 0, "residu": [0, 20], "funciton": 0, "make": [0, 5, 9, 10, 14, 20, 21], "rb": [0, 20, 21], "galerkin_help": 0, "total_basi": 0, "u": [0, 5, 20, 21], "vt": [0, 20, 21], "svd": [0, 5, 20, 21], "full_matric": [0, 20, 21], "fals": [0, 20, 21], "kinetic_ham": 0, "calcul": [0, 9, 10, 20], "strenght": 0, "multipli": [0, 21], "new": [0, 5, 10, 20, 21], "pot_ham": 0, "phi_basi": 0, "l": [0, 10, 20, 21], "reduced_k_ham": 0, "zero": [0, 5, 20, 21], "dtype": 0, "complex_": 0, "reduced_pot_ham": 0, "a0": 0, "three": 0, "princip": [0, 20, 21], "singl": [0, 12], "made": 0, "rbm_stuff": 0, "select": [0, 5, 20, 21], "obtain": [0, 20, 21], "gorgeou": [0, 21], "nbasi": [0, 20, 21], "subplot": [0, 5, 14, 20, 21], "100": [0, 9, 14, 20, 21], "patch": [0, 20, 21], "set_facecolor": [0, 20, 21], "white": [0, 20, 21], "gca": 0, "set_prop_cycl": 0, "set_xlabel": [0, 20, 21], "set_ylabel": [0, 5, 20, 21], "u_": [0, 20, 21], "rm": [0, 20, 21], "train": [0, 9, 10, 20], "solve_linear_system_continuo": 0, "omega_v": 0, "A": [0, 10, 20, 21], "y0": 0, "norm": [0, 5, 20], "sol_rbm": 0, "251": 0, "ms": 0, "157": 0, "\u00b5s": 0, "super": 0, "fast": [0, 10], "specif": 0, "case": [0, 20, 21], "also": [0, 5, 9, 10, 13, 14, 20, 21], "check": [0, 5, 10, 14, 20], "entir": [0, 21], "wf_rbm": 0, "rho_rbm": 0, "rb_rho_movi": 0, "rb_wf_movi": 0, "And": [0, 21], "rb_solver": 0, "877340278": 0, "two": [0, 9, 10, 12, 13, 20], "reproduc": [0, 14], "movie_maker_dens": 0, "list_of_rhos1": 0, "list_of_rhos2": 0, "min": [0, 20], "_densiti": 0, "hf0": 0, "hf": 0, "good": [0, 10, 20], "defin": [0, 5, 13, 20, 21], "metric": 0, "judg": [0, 21], "work": [0, 5, 9, 10, 14, 21], "metric_funct": 0, "true_list_ful": 0, "calculated_list_ful": 0, "x_calc": 0, "x_true": 0, "f_true": 0, "f_calc": 0, "017484555884774013": 0, "For": [0, 5, 9, 12, 14, 20], "had": 0, "befor": [0, 9, 20], "seem": [0, 5], "central": [0, 10], "configur": 0, "order": [0, 5, 9, 10], "cent_config": 0, "truth": 0, "against": 0, "super_hf_config": 0, "20000": 0, "400": 0, "hf_results_sup": 0, "14min": 0, "7s": 0, "34min": 0, "35": 0, "48min": 0, "42": 0, "2min": 0, "26": 0, "It": [0, 21], "veri": [0, 9, 20, 21], "resolut": 0, "like": [0, 1, 9, 13, 20, 21], "hf_configs1": 0, "40": [0, 10], "80": 0, "300": 0, "hf_res_conf1": 0, "comparison": [0, 5], "realli": 0, "precis": 0, "one": [0, 9, 10, 12, 20, 21], "st": 0, "res_dummi": 0, "et": 0, "print": [0, 5, 13, 21], "049285411834717": 0, "19251394271850586": 0, "22269678115844727": 0, "1031908988952637": 0, "5019428730010986": 0, "92787480354309": 0, "22": 0, "826255559921265": 0, "funtion": 0, "evalu": [0, 5, 9, 20], "better": [0, 20], "idea": [0, 1, 20], "rbm_res_mak": 0, "rbm_config": 0, "numb": 0, "numb_basi": 0, "rbm_stuff_temp": 0, "number_calc": 0, "execution_tim": 0, "lambda": [0, 5, 20, 21], "note": [0, 5, 20, 21], "chang": [0, 20, 21], "onli": [0, 9, 20], "expans": [0, 21], "rbm_configs1": 0, "8": [0, 21], "9": [0, 5], "18": [0, 20], "50": [0, 20, 21], "70": 0, "res_rbm1": 0, "res_rb": 0, "002041093859988905": 0, "0016046550153267945": 0, "00186289362000025": 0, "00035783405861791864": 0, "002309614280002279": 0, "0002148565946228041": 0, "0021074782799951207": 0, "295439393898162e": 0, "05": [0, 5], "0021446347699929902": 0, "447686995178741e": 0, "06": [0, 5], "0021209885200005373": 0, "674995713984094e": 0, "003636440180016507": 0, "4941366387066083e": 0, "07": 0, "003813176410003507": 0, "026078519192635e": 0, "08": 0, "0054489675300101225": 0, "585840713840379e": 0, "09": 0, "001217822039998282": 0, "24715577041663786": 0, "0011377165500016417": 0, "021654829164243102": 0, "0014645477900012338": 0, "0021966204719722974": 0, "0017716050700073537": 0, "001162328400873993": 0, "0017636181699890584": 0, "0003315603863323719": 0, "total_rbm": 0, "figsiz": [0, 5, 14], "scatter": [0, 12, 20], "set_xscal": 0, "log": [0, 20, 21], "set_yscal": [0, 20, 21], "23": 0, "rel": 0, "rc": 0, "xtick": 0, "labels": 0, "ytick": 0, "11": 0, "vs": 0, "cat": 0, "excel": [0, 20], "wai": [0, 1, 10, 20, 21], "perform": [0, 9], "other": [0, 9, 10, 13, 20], "each": [0, 5, 9, 20, 21], "trade": 0, "speed": [0, 9, 10], "control": 0, "detail": [0, 13, 21], "mesh": 0, "size": [0, 20, 21], "under": [0, 10], "categor": 0, "observ": [0, 9], "variat": 0, "straighforward": 0, "inde": [0, 9, 20, 21], "give": [0, 5, 9, 10, 21], "your": [0, 1, 5, 9, 13, 14, 20], "own": [0, 20], "model": [0, 5, 10, 21], "close": 0, "might": [0, 5, 20, 21], "pattern": [0, 21], "would": [0, 9, 20], "analysi": [0, 20, 21], "those": [0, 10, 20], "could": [0, 9, 20, 21], "even": [0, 9], "further": [0, 9, 10, 20], "join": [0, 9, 10], "complic": 0, "nonlinear": 0, "eric": [1, 5], "flynn": [1, 5], "content": [1, 13, 20], "daniel": [1, 21], "odel": [1, 21], "ruchi": [1, 9], "garg": [1, 9], "beyer": [1, 20], "webmast": 1, "mani": [1, 5, 9, 10, 13, 21], "everyon": 1, "ha": [1, 5, 10, 21], "book": [1, 9, 10, 13, 14, 20, 21], "add": [1, 20], "anyth": 1, "pleas": [1, 9], "contact": 1, "consid": [5, 10, 20], "align": [5, 14], "phi_": [5, 20, 21], "m": [5, 20], "eigenvalu": 5, "omega": 5, "valu": [5, 9, 20, 21], "without": [5, 10, 20], "repeatli": 5, "directli": 5, "alpha_": 5, "sum_": [5, 20], "a_": [5, 20], "exact": [5, 20, 21], "chosen": [5, 20], "goal": [5, 10], "tutori": [5, 12], "significantli": 5, "cheaper": 5, "interest": [5, 9, 20, 21], "quantifi": 5, "sci": 5, "optim": [5, 20], "special": 5, "demo": 5, "getpsi_x": 5, "definit": [5, 20], "ho": 5, "wavefunct": [5, 20], "zettili": 5, "page": [5, 13], "240": 5, "type": [5, 9], "principl": [5, 9, 20], "se": 5, "mass": [5, 20, 21], "mu": [5, 20, 21], "wf": 5, "1d": 5, "posit": 5, "herm": 5, "hermit": 5, "factori": 5, "25": [5, 20], "getexactlambda": 5, "f": [5, 9, 10, 20, 21], "float": 5, "integ": 5, "frequenc": 5, "2me": 5, "nd": 5, "length": 5, "ndarrai": 5, "construct_h": 5, "2nd": 5, "scheme": [5, 20], "discret": 5, "differenti": [5, 10, 20, 21], "fix": [5, 14, 20, 21], "descript": [5, 9, 10], "2d": 5, "dim": 5, "ident": [5, 20], "toeplitz": 5, "domain": 5, "space": [5, 9, 20], "evect": 5, "eigenvector": 5, "ascend": 5, "eigh": 5, "enumer": 5, "simpson": 5, "getsystem": 5, "psi_arrai": 5, "phi_arrai": 5, "syetem": 5, "projector": [5, 20, 21], "assum": 5, "row": 5, "correspond": 5, "vector": 5, "form": [5, 10, 20, 21], "lambda_": 5, "output": [5, 13], "system": [5, 9, 10, 20, 21], "a_vec": 5, "normal": [5, 20, 21], "fsolv": 5, "choic": [5, 20], "method": [5, 9, 10, 21], "arang": [5, 20, 21], "matmul": 5, "kp": 5, "global": 5, "variabl": [5, 20, 21], "warn": 5, "around": [5, 21], "get": [5, 9, 13, 20, 21], "slow": 5, "x_a": 5, "x_b": 5, "x_arrai": 5, "2001": 5, "These": [5, 9, 20, 21], "reduc": [5, 9, 10, 20], "alpha_v": 5, "hold": [5, 20], "t_eval": 5, "alpha_sampl": 5, "given": [5, 9, 20], "assign": 5, "nth": 5, "sure": [5, 10, 14, 20], "accur": 5, "wf_val": 5, "diff": 5, "ab": [5, 21], "str": 5, "semilog": [5, 20], "set_titl": 5, "sol": [5, 20, 21], "tight_layout": 5, "lambda_exact": 5, "707103656172208": 5, "7071067811865476": 5, "236036727063779": 5, "23606797749979": 5, "1622151589324003": 5, "1622776601683795": 5, "872889593938531": 5, "872983346207417": 5, "1250143395222807e": 5, "125043601093225e": 5, "250123597917323e": 5, "375226888597155e": 5, "appli": [5, 20, 21], "vh": 5, "n_comp": 5, "column": 5, "4": [5, 9, 13, 14, 20, 21], "o": [5, 20], "singular": [5, 20, 21], "sigma": [5, 20], "decreas": [5, 20], "move": [5, 9], "down": 5, "diagon": 5, "decomposit": [5, 20, 21], "tell": 5, "increas": [5, 10], "gain": [5, 9], "less": [5, 20], "inform": [5, 13, 14, 20, 21], "impli": 5, "just": [5, 9, 10, 20, 21], "few": [5, 10, 20], "panel": 5, "know": [5, 20], "rbm": [5, 10, 12, 21], "explicitli": 5, "linearli": [5, 20], "independ": [5, 20, 21], "psi_": 5, "interpret": 5, "enforc": [5, 20], "orthogon": [5, 21], "subspac": 5, "span": 5, "psi": [5, 20, 21], "abil": 5, "arbitrari": 5, "coeffeci": 5, "howev": [5, 20], "unknown": [5, 20], "addit": [5, 9], "come": [5, 10, 21], "alpha_k": 5, "h_k": 5, "same": [5, 9, 20], "abov": [5, 20, 21], "ak_system": 5, "depend": [5, 9, 20, 21], "guess": 5, "excit": [5, 12], "mai": 5, "rest": [5, 13, 20], "ok": 5, "ak_guess": 5, "a_sol": 5, "appproxim": 5, "approx_sol": 5, "true_wf": 5, "02001628": 5, "72240374": 5, "01421624": 5, "31046244": 5, "001175561288746": 5, "0011755612887460742": 5, "edit": 9, "dive": 9, "basi": [9, 10, 20], "must": [9, 20], "inward": 9, "deep": 9, "insid": 9, "soul": 9, "ask": 9, "fundament": 9, "question": 9, "learn": 9, "practition": 9, "view": 9, "thought": 9, "algorithm": [9, 20], "mimic": 9, "expens": 9, "littl": 9, "loss": 9, "gigant": 9, "too": [9, 10], "hassl": 9, "should": [9, 13, 20], "alwai": [9, 21], "neglig": 9, "At": 9, "dai": 9, "time": [9, 20, 21], "faster": 9, "great": [9, 20], "win": 9, "almost": 9, "everi": [9, 10, 20], "research": [9, 10], "group": 9, "walk": 9, "short": [9, 20], "benefici": 9, "decis": 9, "amount": 9, "scienc": [9, 10, 21], "imagin": 9, "colleagu": 9, "develop": [9, 10], "nice": [9, 20], "list": [9, 12, 20], "explain": [9, 20, 21], "phenomena": 9, "natur": [9, 10], "neutron": [9, 21], "nucleu": [9, 20], "beam": 9, "strongli": 9, "figure1": 9, "hypothet": 9, "situat": [9, 21], "yourself": 9, "read": [9, 10, 20], "red": 9, "compar": 9, "blue": 9, "friend": [9, 10], "laboratori": 9, "measur": 9, "campaign": 9, "respons": 9, "until": 9, "match": 9, "perhap": [9, 10], "hundr": 9, "thousand": 9, "task": 9, "becom": [9, 10], "consum": 9, "hypotet": 9, "opposit": 9, "side": 9, "burden": [9, 10], "spectrum": 9, "instead": [9, 20], "bayesian": 9, "posterior": 9, "distribut": [9, 20], "thing": [9, 13, 20], "quickli": 9, "out": [9, 14, 20, 21], "hand": [9, 20], "band": 9, "wa": [9, 10, 21], "built": [9, 13], "sampl": [9, 14, 20], "000": 9, "sometim": 9, "million": 9, "were": [9, 20, 21], "parallel": 9, "cluster": 9, "turn": 9, "year": [9, 10], "simpl": [9, 10, 12, 20], "millennia": 9, "requir": 9, "resourc": [9, 10], "quantif": [9, 10], "within": 9, "lifetim": 9, "put": [9, 20, 21], "hero": 9, "carefulli": 9, "glori": 9, "consist": 9, "synchron": 9, "perfect": 9, "provid": [9, 20], "reliabl": 9, "scientif": [9, 10], "commun": [9, 10], "facil": 9, "investig": 9, "rare": 9, "nuclei": [9, 21], "secar": 9, "separ": [9, 20], "captur": 9, "reaction": 9, "dedic": 9, "cross": [9, 20], "probabl": [9, 21], "nuclear": [9, 10], "product": [9, 20, 21], "unreact": 9, "seri": [9, 10], "magnet": 9, "illustr": [9, 20], "ion": [9, 14], "path": [9, 13], "thei": [9, 21], "travers": 9, "between": [9, 20, 21], "target": [9, 20], "locat": [9, 20], "detect": 9, "plane": 9, "adapt": 9, "focus": 9, "defocus": 9, "rai": 9, "pass": [9, 20], "dipol": 9, "quadrupol": 9, "hexapol": 9, "octupol": 9, "wien": 9, "filter": 9, "steer": 9, "desir": [9, 20], "manner": 9, "challeng": [9, 10], "involv": [9, 20], "computation": 9, "simul": 9, "manual": 9, "fine": 9, "effect": [9, 20], "larg": [9, 20], "help": 9, "tremend": 9, "process": 9, "reduct": [9, 10, 20], "greatli": 9, "dure": [9, 10], "valuabl": 9, "therebi": 9, "spent": 9, "apparatu": 9, "maxim": [9, 20], "activ": 9, "collect": 9, "hope": [9, 10], "serv": [9, 10, 20], "appet": 9, "motiv": [9, 20], "techniqu": [9, 10, 20], "intuit": [9, 20], "nutshel": 9, "vari": [9, 20, 21], "highli": 9, "redund": 9, "As": [10, 14, 20, 21], "scientist": 10, "degre": 10, "freedom": 10, "theori": 10, "includ": [10, 13, 14], "ever": 10, "grow": 10, "been": 10, "art": 10, "than": [10, 20, 21], "featur": 10, "success": 10, "share": 10, "centuri": 10, "novel": 10, "approach": [10, 21], "mostli": 10, "decad": 10, "power": [10, 20], "area": 10, "allow": [10, 21], "jupyt": [10, 13, 14], "overview": 10, "umbrella": 10, "introduc": 10, "ago": 10, "constantli": [10, 21], "improv": 10, "especi": 10, "There": [10, 14], "aricl": 10, "subject": 10, "unfortun": 10, "literatur": [10, 20], "strong": 10, "background": 10, "formal": 10, "abstract": 10, "mathemat": 10, "treatment": 10, "concept": 10, "weak": 10, "formul": 10, "kolmogorov": 10, "width": 10, "coercit": 10, "continu": [10, 20], "discuss": 10, "decid": 10, "audienc": 10, "broad": 10, "physic": 10, "lean": 10, "bit": 10, "toward": 10, "undergradu": 10, "graduat": 10, "student": 10, "matur": 10, "who": 10, "tool": 10, "favor": 10, "concret": 10, "while": [10, 20], "connect": 10, "pursu": 10, "rout": 10, "extens": [10, 20, 21], "dirac": 10, "bra": 10, "ket": 10, "common": 10, "familiar": 10, "accompani": 10, "explan": 10, "code": [10, 12, 13, 20, 21], "written": [10, 13, 21], "plai": 10, "fun": 10, "enjoi": 10, "topic": [10, 20], "implement": [10, 20, 21], "support": [10, 13], "amaz": 10, "discoveri": 10, "ongo": 10, "progress": 10, "updat": 10, "email": 10, "suggest": 10, "discov": 10, "bug": 10, "recur": 10, "journal": 10, "club": 10, "team": 10, "wonder": 10, "driven": 10, "engin": 10, "video": 10, "brunton": 10, "kutz": 10, "quarteroni": 10, "manzoni": 10, "negri": 10, "certifi": 10, "parametr": 10, "hesthaven": 10, "rozza": 10, "stamm": 10, "survei": 10, "base": [10, 13], "p": [10, 20, 21], "benner": 10, "gugercin": 10, "willcox": 10, "melendez": 10, "c": [10, 21], "drischler": 10, "furnstahl": 10, "garcia": 10, "zhang": 10, "soon": [10, 12], "articl": 10, "bodi": [10, 12], "lee": 10, "bay": 10, "goe": [10, 21], "uncertainti": 10, "covari": 10, "vien": 10, "piekarewicz": 10, "fill": 12, "applic": 12, "paper": [12, 20, 21], "kylegodbei": 12, "contributor": 12, "found": 12, "channel": 12, "empir": [12, 21], "affin": 12, "decompos": 12, "document": [13, 14], "instruct": 13, "off": 13, "With": 13, "direct": [13, 20], "block": 13, "default": 13, "kernel": 13, "displai": 13, "jupytext": 13, "convert": 13, "understand": 13, "top": 13, "presenc": 13, "That": 13, "treat": 13, "markdownfil": 13, "md": 13, "emb": 14, "html": 14, "etc": 14, "post": 14, "add_": 14, "mbox": 14, "la_": 14, "tex": 14, "But": [14, 20], "escap": 14, "dollar": 14, "sign": 14, "keep": [14, 20], "guid": 14, "rcparam": 14, "cycler": 14, "_ioncontext": 14, "0x7f9e04f21290": 14, "random": [14, 20, 21], "seed": 14, "19680801": 14, "logspac": 14, "randn": 14, "ii": 14, "cmap": 14, "cm": 14, "coolwarm": 14, "prop_cycl": 14, "line2d": 14, "custom_lin": 14, "cold": 14, "medium": 14, "hot": 14, "inteprol": 20, "known": 20, "eim": 20, "radial": [20, 21], "dr": [20, 21], "ell": [20, 21], "wood": 20, "saxon": 20, "v_0": [20, 21], "correlti": 20, "momentum": [20, 21], "scale": [20, 21], "tild": [20, 21], "z": 20, "ve": 20, "notic": 20, "distinc": 20, "high": 20, "fidel": 20, "vario": 20, "alpha_i": 20, "v_": [20, 21], "s_i": [20, 21], "z_i": 20, "p_i": 20, "alpha_t": 20, "re": [20, 21], "construct": [20, 21], "opper": 20, "clever": 20, "behind": 20, "beta_k": 20, "u_k": 20, "determin": 20, "somehow": 20, "radiu": 20, "ommit": 20, "forward": 20, "analog": 20, "did": [20, 21], "equiv": 20, "free": [20, 21], "suitabl": 20, "encapsul": 20, "associ": 20, "becaus": [20, 21], "don": 20, "access": 20, "incur": 20, "cost": 20, "therefor": 20, "beta": 20, "exactli": 20, "equal": [20, 21], "origin": [20, 21], "s_j": 20, "repres": 20, "expect": [20, 21], "reason": [20, 21], "constantn": 20, "hbarc": [20, 21], "197": [20, 21], "mev": [20, 21], "fm": [20, 21], "939": [20, 21], "center": [20, 21], "woods_saxon": 20, "v0": [20, 21], "woods_saxon_tild": 20, "test": [20, 21], "evenli": 20, "inexpens": 20, "simpli": 20, "r_min": 20, "r_max": 20, "generate_rand_log": 20, "utrain": 20, "spit": 20, "outer": 20, "shape": [20, 21], "tt": 20, "ss": 20, "s_endpt": [20, 21], "1e": [20, 21], "dimensionless": [20, 21], "2000": [20, 21], "umat": 20, "ylabel": 20, "healthi": 20, "diffract": 20, "minima": 20, "elast": 20, "quit": 20, "sensit": 20, "word": 20, "pca": 20, "sv": 20, "index": [20, 21], "sing": 20, "val": 20, "ratio": 20, "3f": 20, "format": 20, "pair": 20, "place": 20, "varianc": 20, "highest": 20, "standard": 20, "deviat": 20, "var": 20, "axi": [20, 21], "s_k": 20, "uniformli": 20, "below": 20, "altern": 20, "pick": 20, "downsid": 20, "nois": 20, "averag": 20, "accord": 20, "sequenti": 20, "essenc": 20, "greedi": 20, "cover": 20, "ps": 20, "trapz": 20, "sum": [20, 21], "sample_from_vari": 20, "sample_pt": 20, "els": 20, "96301274": 20, "47062382": 20, "49801361": 20, "37977634": 20, "8075444": 20, "23671835": 20, "onc": 20, "anoth": [20, 21], "simultan": 20, "mathbf": 20, "overlin": 20, "jk": 20, "u_j": 20, "unkown": 20, "interp": 20, "ki": 20, "marker": 20, "einsum": 20, "inv": 20, "rk": 20, "magnitud": 20, "saniti": 20, "ws_affine_decomp": 20, "uad": 20, "r_test": 20, "By": 20, "mess": 20, "extrapol": 20, "r_": 20, "notin": 20, "max": [20, 21], "visual": 20, "poor": 20, "agreement": 20, "displac": 20, "extrem": 20, "infinit": 20, "smaller": 20, "borrow": 20, "machineri": 20, "kei": 20, "initial_condit": [20, 21], "solve_s": [20, 21], "theta": [20, 21], "rtol": [20, 21], "atol": [20, 21], "dense_output": [20, 21], "n_train_pt": 20, "training_set": 20, "training_point": [20, 21], "r_train": 20, "training_solut": [20, 21], "soln": 20, "subtract": 20, "avoid": 20, "divid": 20, "exponenti": [20, 21], "training_solutions_sub": [20, 21], "zeros_lik": [20, 21], "preliminari": 20, "u_rb": 20, "s_rb": 20, "0x7efef7f99a90": 20, "phi_i": 20, "itself": 20, "psi_i": 20, "psi_j": [20, 21], "ds": [20, 21], "inhomogen": 20, "vec": 20, "act": [20, 21], "psi_1": 20, "psi_2": 20, "psi_n": 20, "a_2": [20, 21], "kappa": 20, "beta_": 20, "Then": 20, "grab": 20, "weight": 20, "b_j": 20, "rh": 20, "lh": 20, "pre": 20, "someth": 20, "u_proj": 20, "si": 20, "sj": 20, "sk": 20, "u_0": 20, "precomput": 20, "leav": 20, "angular": [20, 21], "gradient": 20, "minu": 20, "care": 20, "rbm_emul": 20, "appropri": 20, "onto": [20, 21], "a_pot": 20, "b_pot": 20, "emu_proj": 20, "9454773": 20, "9964038": 20, "32121418": 20, "39575502": 20, "09115662": 20, "back": 20, "emu": [20, 21], "Not": 20, "bad": 20, "asymptot": 21, "optic": 21, "introduct": 21, "pr": 21, "possibl": 21, "proce": 21, "minnesota": 21, "0r": 21, "0s": 21, "487r": 21, "465r": 21, "kept": 21, "spars": 21, "style": 21, "coupl": 21, "aesthet": 21, "gr": 21, "v0r": 21, "91": 21, "85": 21, "mn_potenti": 21, "v_0r": 21, "487": 21, "465": 21, "mn_potential_tild": 21, "rho": 21, "kr": 21, "ns": 21, "s_mesh": 21, "119": 21, "51219512195122": 21, "14": 21, "634146341463415": 21, "139": 21, "02439024390245": 21, "878048780487805": 21, "158": 21, "53658536585365": 21, "48": 21, "78048780487805": 21, "178": 21, "0487804878049": 21, "117": 21, "07317073170732": 21, "5609756097561": 21, "131": 21, "70731707317074": 21, "217": 21, "0731707317073": 21, "126": 21, "82926829268293": 21, "236": 21, "58536585365854": 21, "82": 21, "92682926829268": 21, "256": 21, "0975609756098": 21, "175": 21, "609756097561": 21, "275": 21, "19": 21, "295": 21, "1219512195122": 21, "170": 21, "73170731707316": 21, "decai": 21, "indic": 21, "succes": 21, "higher": 21, "outsid": 21, "region": 21, "express": 21, "combin": 21, "relat": 21, "sine": 21, "cosin": 21, "why": 21, "particularli": 21, "overal": 21, "zig": 21, "zag": 21, "stair": 21, "four": 21, "interestingli": 21, "homogen": 21, "valid": 21, "impos": 21, "similar": 21, "avaialbl": 21, "abscenc": 21, "utild": 21, "newaxi": 21, "toarrai": 21, "carri": 21, "aa": 21, "inner": 21, "psi_k": 21, "sens": 21, "sourc": 21, "Near": 21, "endpoint": 21, "signific": 21, "a_right": 21, "55324541e": 21, "00": 21, "69173226e": 21, "01": 21, "70915793e": 21, "42882188e": 21, "03": 21, "enough": 21, "biggest": 21, "exmapl": 21, "think": 21, "doesn": 21}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"applic": [0, 2, 3, 5, 6, 20, 21], "4": [0, 2], "time": 0, "depend": 0, "system": 0, "evolut": 0, "reduc": [0, 12, 15, 18, 21], "space": 0, "high": [0, 21], "fidel": [0, 21], "solver": [0, 21], "approach": [0, 7, 8], "basi": [0, 5, 12, 15, 18, 19, 21, 22], "method": [0, 12, 18, 20], "compar": [0, 5], "rbm": [0, 7, 20], "fix": 0, "valu": 0, "omega": 0, "contributor": 1, "nuclear": [2, 12, 21], "dft": 2, "2": [3, 21], "The": [3, 4, 5, 6, 7, 8, 18, 20], "gross": 3, "pitaevskii": 3, "equat": 3, "greedi": 4, "algorithm": 4, "1": [5, 6], "quantum": [5, 6], "harmon": [5, 6], "oscil": [5, 6], "background": 5, "solv": 5, "exactli": 5, "visit": 5, "pca": [5, 16], "doctor": 5, "construct": 5, "ritz": 8, "variat": 8, "why": 9, "emul": 9, "first": 9, "case": 9, "calibr": 9, "quantifi": 9, "uncertainti": 9, "physic": [9, 12], "model": 9, "second": 9, "experiment": 9, "optim": 9, "onlin": 9, "tune": 9, "up": 9, "introduct": [10, 12], "redund": 11, "simul": 11, "notebook": [13, 14], "myst": [13, 14], "markdown": [13, 14], "an": 13, "exampl": 13, "cell": 13, "creat": 13, "quickli": 13, "add": 13, "yaml": 13, "metadata": 13, "content": 14, "code": 14, "block": 14, "output": 14, "lagrang": 15, "pod": 16, "find": 17, "coeffici": 17, "project": 17, "choos": 19, "train": [19, 21], "3": 20, "empir": 20, "interpol": 20, "build": 20, "affin": 20, "decompos": 20, "problem": 20, "two": 21, "bodi": 21, "singl": 21, "channel": 21, "scatter": 21, "modifi": 22, "stretch": 22, "translat": 22}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})