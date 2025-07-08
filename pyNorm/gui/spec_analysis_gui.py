import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib as mpl
import pickle
import glob
import os
import sys
from PyQt5.QtWidgets import QMessageBox, QApplication, QDialog, QVBoxLayout, QLabel, QPushButton, QSizePolicy, QScrollArea, QWidget
from PyQt5.QtWidgets import QTableWidgetItem
import PyQt5.QtWidgets as qtw
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
from astropy.io import ascii
from astropy.table import Table
from pyND.absorption import logmean  # Assuming this is your custom module
from concurrent.futures import ThreadPoolExecutor
from functools import partial  # For optimized callbacks
import traceback
import matplotlib as mpl
from cycler import cycler
NAV_COLORS = [
    '#1f77b4',  # dark blue
    '#ff7f0e',  # dark orange
    '#2ca02c',  # dark green
    '#6816B4',  # dark purple
    '#5A3D38',  # brown
    '#B4007D',  # dark pink
    '#7f7f7f',  # dark gray
    '#0095A4',  # cyan
    '#1900FF',  # olive
]
mpl.rcParams['axes.prop_cycle'] = cycler('color', NAV_COLORS)
def f_lambda(wavelength, f_value):
    return wavelength * f_value

def average_log_errors(logN, err_lo, err_hi):
    """Convert asymmetric log errors to a symmetric error in linear space, then back to log space."""
    linear = 10**logN
    sig_lo = linear - 10**(logN - err_lo)
    sig_hi = 10**(logN + err_hi) - linear
    return np.mean([sig_lo, sig_hi])

def correct_saturation(logN, err_logN, corr_max = 0.15, printresults=True):
    """Correct mildly-saturated Na(v) results for doublets following the recommendation in Savage & Sembach (1991). **This assumes a factor of 2x difference in f-values.

    CALLING SEQUENCE:
    sscorrect, logn_array, err(logn)_array

    INPUTS:
    logN     -- a two-element array holding the Na(v) values [strong, weak].
    err_logN -- a two-element array holding the Na(v) errors [strong, weak].

    OPTIONAL INPUTS:
    corr_max = 0.15     -- The assumed maximum value of valid corrections.
    printresults = True -- print results to the console.

    OUTPUTS:
    logNf     -- Final corrected column density.
    err_logNf -- Error in final column density (-2 == saturated).
    """

    # The difference in columns is
    diff=logN[1]-logN[0]
    # Check that the column densities aren't reversed.
    # Assume they are if the difference in columns is negative.
    if diff < 0:
        print, "The logN difference is negative! Did you put the strong line first? If so, there may be a problem, BUT ASSUMING THERE IS NONE...."

        logN = logN[::-1]
        err_logN = err_logN[::-1]
        diff=logN[1]-logN[0]

    # Calculate the correction using polynomial fit to SS1991 results.
    ssdiff=diff
    sscorrect=16.026*ssdiff**3 - 0.507*ssdiff**2 \
                + 0.9971*ssdiff + 5.e-5


    err_logNf = np.sqrt(err_logN[1]**2 + \
        ((48.078*ssdiff**2 - 1.014*ssdiff + 0.9971)*np.sqrt(err_logN[1]**2 \
        + err_logN[0]**2) )**2 )

    if sscorrect >= corr_max:
        print("Your correction exceeded the maximum correction, d(logN)_max = {0:0.3f}.".format(corr_max))
        print("Applying the maximum correction and assuming the final result is saturated.")

        logNf = logN[1]+corr_max
        err_logNf = -2
        if printresults:
            print("Final result: logN_final > {0:0.3f}".format(logNf))
    else:
        logNf=logN[1]+sscorrect
        if printresults:
            print("Final result: logN_final = {0:0.3f}+/-{1:0.3f}".format(logNf,err_logNf))

    return logNf, err_logNf

class SpectroscopicAnalysisGUI(qtw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Spectroscopic Absorption Line Analysis")
        self.setGeometry(100, 100, 1800, 1000)
        
        # Data structures
        self.target_data = {}
        self.spec_data = {}
        self.summary = None
        self.va_label = None
        self.va_err_label = None
        self.current_ion = ""
        self.current_target = ""
        self.current_redshift = 0.0
        self.ion_list = []
        self.ions = []
        self.line_vars = {}  # For line checkboxes
        self.quality_vars = {}
        self.entry_widgets = {}

        # Set matplotlib parameters
        mpl.rcParams['mathtext.default'] = 'it'
        mpl.rcParams['font.size'] = 12
        mpl.rc('font', family='Times New Roman')
        mpl.rc('xtick', labelsize=14)
        mpl.rc('ytick', labelsize=14)
        mpl.rc('axes', labelsize=14)
        
        # Create GUI
        self.setup_ui()
        
        # Check for existing summary file
        self.check_summary_file()
        
        # Cache for plot data to improve performance
        self.plot_cache = {}
    
    def setup_ui(self):
        # Central widget
        central_widget = qtw.QWidget()
        self.setCentralWidget(central_widget)
        main_layout = qtw.QVBoxLayout(central_widget)
        
        # Top control bar
        top_frame = qtw.QHBoxLayout()
        main_layout.addLayout(top_frame)
        
        # Load button
        self.load_button = qtw.QPushButton("Load Data Files")
        self.load_button.clicked.connect(self.load_data_files)
        top_frame.addWidget(self.load_button)
        
        # Target selection
        top_frame.addWidget(qtw.QLabel("Target:"))
        self.target_dropdown = qtw.QComboBox()
        self.target_dropdown.currentTextChanged.connect(self.target_changed)
        top_frame.addWidget(self.target_dropdown)
        
        # Ion selection
        top_frame.addWidget(qtw.QLabel("Ion:"))
        self.ion_dropdown = qtw.QComboBox()
        self.ion_dropdown.setMaximumWidth(100)
        self.ion_dropdown.currentTextChanged.connect(self.ion_changed)
        top_frame.addWidget(self.ion_dropdown)
        
        top_frame.addStretch()
        
        # Save button
        self.save_button = qtw.QPushButton("Save Results")
        self.save_button.clicked.connect(self.save_results)
        self.save_button.setEnabled(False)
        top_frame.addWidget(self.save_button)
        
        self.stacked_plots_button = qtw.QPushButton("Stacked Plots")
        self.stacked_plots_button.clicked.connect(self.show_stacked_plots)
        self.stacked_plots_button.setEnabled(False)
        top_frame.addWidget(self.stacked_plots_button)

        # Create splitter (equivalent to PanedWindow)
        splitter = qtw.QSplitter(qtc.Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left side - Plot area
        plot_widget = qtw.QWidget()
        plot_layout = qtw.QVBoxLayout(plot_widget)
        
        # Create matplotlib figure and canvas
        self.figure = plt.Figure(figsize=(8, 6), dpi=100)
        self.ax = self.figure.add_subplot(111)
        
        self.canvas = FigureCanvas(self.figure)
        plot_layout.addWidget(self.canvas)
        
        # Add toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        
        splitter.addWidget(plot_widget)
        
        # Right side - Controls and data
        controls_widget = qtw.QWidget()
        controls_layout = qtw.QVBoxLayout(controls_widget)
        
        # Line selection group
        font = qtg.QFont("Times New Roman", 16)
        
        lines_group = qtw.QGroupBox("Available Absorption Lines")
        lines_group.setFont(font)
        lines_group.setFixedHeight(200)  # scroll height + padding
        lines_layout = qtw.QVBoxLayout(lines_group)
        
        # Create scroll area for checkboxes
        self.lines_scroll = qtw.QScrollArea()
        self.lines_scroll.setMaximumHeight(180)  # or adjust based on how many rows you expect
        self.lines_scroll.setWidgetResizable(True)
        self.lines_scroll_widget = qtw.QWidget()
        self.lines_scroll_layout = qtw.QVBoxLayout(self.lines_scroll_widget)
        self.lines_scroll.setWidget(self.lines_scroll_widget)
        lines_layout.addWidget(self.lines_scroll)
        
        controls_layout.addWidget(lines_group)

        # Parameters group
        params_group = qtw.QGroupBox("Measurement Parameters")
        params_group.setFont(font)
        params_layout = qtw.QGridLayout(params_group)
    
        # Column density display
        params_layout.addWidget(qtw.QLabel("Column Density (N):"), 0, 0)
        self.n_value_label = qtw.QLabel("nan")
        params_layout.addWidget(self.n_value_label, 0, 1)
        
        params_layout.addWidget(qtw.QLabel("Error (low):"), 1, 0)
        self.n_err_low_label = qtw.QLabel("nan")
        params_layout.addWidget(self.n_err_low_label, 1, 1)
        
        params_layout.addWidget(qtw.QLabel("Error (high):"), 2, 0)
        self.n_err_high_label = qtw.QLabel("nan")
        params_layout.addWidget(self.n_err_high_label, 2, 1)
        
        # Reliability flag
        params_layout.addWidget(qtw.QLabel("Reliability Flag:"), 3, 0)
        self.reliability_combo = qtw.QComboBox()
        self.reliability_combo.addItems([
            "1 - Non-detection / Multiple lines", 
            "2 - One detection consistent with UL", 
            "3 - Single line / Several saturated"
        ])
        params_layout.addWidget(self.reliability_combo, 3, 1)
        params_layout.addWidget(qtw.QLabel("Redshift Error:"), 4, 0)
        self.z_err_label = qtw.QLabel("nan")
        params_layout.addWidget(self.z_err_label, 4, 1)

        # Detection flag
        params_layout.addWidget(qtw.QLabel("Detection Flag:"), 5, 0)
        self.detection_combo = qtw.QComboBox()
        self.detection_combo.addItems([
            "0 - Detection", 
            "-1 - Upper limit/Non-detection", 
            "-2 - Lower limit/Saturation"
        ])
        params_layout.addWidget(self.detection_combo, 5, 1)
        
        controls_layout.addWidget(params_group)
        
        # Action buttons
        button_layout = qtw.QHBoxLayout()
        
        self.update_button = qtw.QPushButton("Update Measurement")
        self.update_button.clicked.connect(self.update_measurement)
        self.update_button.setEnabled(False)
        button_layout.addWidget(self.update_button)
        
        self.next_button = qtw.QPushButton("Next Ion")
        self.next_button.clicked.connect(self.next_ion)
        self.next_button.setEnabled(False)
        button_layout.addWidget(self.next_button)
        
        controls_layout.addLayout(button_layout)
        
        # Summary Table
        table_group = qtw.QGroupBox("Summary Table")
        table_group.setFont(font)

        table_layout = qtw.QVBoxLayout(table_group)
        
        # Create table widget
        self.summary_table = qtw.QTableWidget()
        self.summary_table.setColumnCount(11)
        self.summary_table.setHorizontalHeaderLabels([
            "Redshift", "Ion", "va", "va_err", "z_err", "N", "N_sig_lo", "N_sig_hi", 
            "Corrected?", "Reliability", "Detection"
                ])

        
        # Set column widths
        header = self.summary_table.horizontalHeader()
        header.setStretchLastSection(True)
        self.summary_table.setColumnWidth(0, 100)
        self.summary_table.setColumnWidth(1, 100)
        self.summary_table.setColumnWidth(2, 80)
        self.summary_table.setColumnWidth(3, 100)
        self.summary_table.setColumnWidth(4, 100)
        self.summary_table.setColumnWidth(5, 100)
        self.summary_table.setColumnWidth(6, 100)
        self.summary_table.setColumnWidth(7, 80)
        self.summary_table.setColumnWidth(8, 80)
        
        table_layout.addWidget(self.summary_table)
        controls_layout.addWidget(table_group)
        
        splitter.addWidget(controls_widget)
        
        # Set splitter proportions
        splitter.setSizes([900, 900])

    def show_stacked_plots(self):
        """Open the stacked plots dialog"""
        if not self.spec_data:
            qtw.QMessageBox.warning(self, "Warning", "No spectral data loaded. Please select a target first.")
            return
    
        self.stacked_plots_dialog = StackedPlotsDialog(self)
        self.stacked_plots_dialog.show()
        self.save_stacked_plots()

    def save_stacked_plots(self):
        """Save all stacked plots to the target folder"""
        if not self.spec_data:
            return

        try:
        # Get target name and create filename-safe version
            target_name = self.current_target
            if not target_name:
                target_name = "unknown_target"
        
        # Remove any characters that might be problematic in filenames
            safe_target_name = "".join(c for c in target_name if c.isalnum() or c in (' ', '-', '_')).rstrip()
            safe_target_name = safe_target_name.replace(' ', '_')
        
        # Get redshift for filename and folder detection
            redshift = getattr(self, 'current_redshift', 'unknown')
        
        # Determine the target folder using the same logic as load_data_files
            target_folder = None
        
        # Method 1: Use current_file_path if available
            if hasattr(self, 'current_file_path') and self.current_file_path:
                target_folder = os.path.dirname(self.current_file_path)
        
        # Method 2: Search in the standard directory structure
            if not target_folder or not os.path.exists(target_folder):
                super_dirs = ["Metals_JMO13", "Metals_R11"]
            
                for super_dir in super_dirs:
                    if os.path.exists(super_dir):
                    # Look for the target directory pattern
                        target_pattern = os.path.join(super_dir, "*", f"results_z_{redshift}")
                        matching_dirs = glob.glob(target_pattern)
                    
                    # Also try with wildcards in case redshift format varies
                        if not matching_dirs:
                            target_pattern = os.path.join(super_dir, f"*{target_name}*", "results_z_*")
                            matching_dirs = glob.glob(target_pattern)
                        
                        if matching_dirs:
                            target_folder = matching_dirs[0]  # Use the first match
                            break
        
        # Method 3: Fallback to current working directory
            if not target_folder or not os.path.exists(target_folder):
                target_folder = os.getcwd()
                print(f"Warning: Could not find target-specific folder, saving to current directory: {target_folder}")
        
        # Create base filename
            base_filename = f"{safe_target_name}_z{redshift}_stacked_plots"
        
        # Save plots for each tab
            saved_files = []
        
        # Save All Transitions plot
            if hasattr(self.stacked_plots_dialog, 'all_transitions_fig'):
                filename = os.path.join(target_folder, f"{base_filename}_all_transitions.pdf")
                self.stacked_plots_dialog.all_transitions_fig.savefig(filename, dpi=300, bbox_inches='tight')
                saved_files.append(filename)
        
        # Save HI Transitions plot
            if hasattr(self.stacked_plots_dialog, 'hi_transitions_fig'):
                filename = os.path.join(target_folder, f"{base_filename}_hi_transitions.pdf")
                self.stacked_plots_dialog.hi_transitions_fig.savefig(filename, dpi=300, bbox_inches='tight')
                saved_files.append(filename)
            
        # Save Doublets plot
            if hasattr(self.stacked_plots_dialog, 'doublets_fig'):
                filename = os.path.join(target_folder, f"{base_filename}_doublets.pdf")
                self.stacked_plots_dialog.doublets_fig.savefig(filename, dpi=300, bbox_inches='tight')
                saved_files.append(filename)
        
            if saved_files:
                print(f"Stacked plots saved to: {target_folder}")
                for file in saved_files:
                    print(f"  - {os.path.basename(file)}")
            else:
                print("No plots were saved (no figures found)")
    
        except Exception as e:
            print(f"Error saving stacked plots: {e}")
            qtw.QMessageBox.warning(self, "Warning", f"Failed to save plots: {str(e)}")
    
    def check_summary_file(self):
        """Check if summary file exists and load it"""
        os.makedirs('Summary', exist_ok=True)
        if os.path.exists('Summary/summary.csv'):
            try:
                self.summary = Table.read('Summary/summary.csv', format='csv')
                self.display_summary_table()
            except Exception as e:
                qtw.QMessageBox.warning(self, "Warning", f"Error loading summary file: {str(e)}")
                self.summary = Table(names=['redshift', 'targname', 'ion', 'va','va_err', 'z_err', 'N', 'N_sig_lo', 'N_sig_hi', 'corrected_flag', 'reliability', 'detection'],
                    dtype=[float, str, str, float, float, float, float, float, float, int, int, int])

        else:
            self.summary = Table(names=['redshift', 'targname', 'ion', 'va','va_err', 'z_err', 'N', 'N_sig_lo', 'N_sig_hi', 'corrected_flag', 'reliability', 'detection'],
                    dtype=[float, str, str, float, float, float, float, float, float, int, int, int])
    
    def load_data_files(self):
        """Load data files from standard directories"""
        # Show busy cursor
        qtw.QApplication.setOverrideCursor(qtc.Qt.WaitCursor)
        
        try:
            # Look under both superdirectories
            super_dirs = ["Metals_JMO13", "Metals_R11"]
            file_list = []

            for super_dir in super_dirs:
                if os.path.exists(super_dir):
                    pattern = os.path.join(super_dir, "*", "results_z_*", "*_Merged.pkl")
                    file_list.extend(glob.glob(pattern))

            if not file_list:
                qtw.QMessageBox.warning(self, "Warning", "No data files found in Metals_JMO13 or Metals_R11")
                return

            # Load the first file to populate the UI
            self.load_targets(file_list)
        
        finally:
            # Restore cursor
            qtw.QApplication.restoreOverrideCursor()
    
    @staticmethod
    def load_pickle_wrapper(file):
        try:
            with open(file, 'rb') as f:
                spec = pickle.load(f)
            ion_list = list(spec.keys())
            if ion_list:
                targname = str(spec['Sightline']).strip()
                redshift = float(spec[ion_list[0]]['z'])
            # Ensure consistent formatting
                key = f"{targname} (z={redshift:.5f})"
                return key, (file, targname, redshift)

        except Exception as e:
            print(f"Error loading {file}: {str(e)}")
        return None

    def load_targets(self, file_list):
        """Load target information from files and populate dropdowns"""
        targets = {}
        with ThreadPoolExecutor(max_workers=8) as executor:
            results = list(executor.map(SpectroscopicAnalysisGUI.load_pickle_wrapper, file_list))
    
        for result in results:
            if result:
                key, value = result
                targets[key] = value
    
        if targets:
        # IMPORTANT: Set target_data BEFORE populating dropdown
            self.target_data = targets
        
        # Temporarily disconnect the signal to prevent premature triggering
            self.target_dropdown.currentTextChanged.disconnect()
        
        # Now populate target dropdown
            target_names = sorted(list(targets.keys()))
            self.target_dropdown.clear()
            self.target_dropdown.addItems(target_names)
        
        # Reconnect the signal
            self.target_dropdown.currentTextChanged.connect(self.target_changed)
        
        # Enable appropriate buttons
            self.save_button.setEnabled(True)
            self.stacked_plots_button.setEnabled(True)
        # Now manually trigger the first selection
            if target_names:
                self.target_dropdown.setCurrentIndex(0)
                self.target_changed()  # Now call it manually when data is ready
        else:
            qtw.QMessageBox.warning(self, "Warning", "No valid data found in selected files")

    def debug_target_keys(self):
        """Debug method to print target keys and dropdown text"""
        print("=== TARGET DATA KEYS ===")
        for i, key in enumerate(self.target_data.keys()):
            print(f"{i}: '{key}' (len={len(key)})")
    
        print("\n=== DROPDOWN ITEMS ===")
        for i in range(self.target_dropdown.count()):
            item_text = self.target_dropdown.itemText(i)
            print(f"{i}: '{item_text}' (len={len(item_text)})")
    
        print(f"\n=== CURRENT SELECTION ===")
        current = self.target_dropdown.currentText()
        print(f"Current: '{current}' (len={len(current)})")

    def target_changed(self):
        qtw.QApplication.setOverrideCursor(qtc.Qt.WaitCursor)

        try:
            dropdown_text = self.target_dropdown.currentText().strip()
        
        # Try exact match first
            matched_key = None
            for key in self.target_data:
                if key.strip() == dropdown_text:
                    matched_key = key
                    break
        
        # If no exact match, try case-insensitive match
            if matched_key is None:
                for key in self.target_data:
                    if key.strip().lower() == dropdown_text.lower():
                        matched_key = key
                        break
        
        # If still no match, try to find partial matches
            if matched_key is None:
                potential_matches = []
                for key in self.target_data:
                    if dropdown_text in key or key.strip() in dropdown_text:
                        potential_matches.append(key)
            
                if len(potential_matches) == 1:
                    matched_key = potential_matches[0]
                elif len(potential_matches) > 1:
                    print(f"Multiple potential matches found: {potential_matches}")
                    matched_key = potential_matches[0]  # Use first match

            if matched_key is None:
                print(f"ERROR: Target '{dropdown_text}' not found in loaded data.")
                print("Available targets:")
                for key in self.target_data.keys():
                    print(f"  '{key}'")
                qtw.QMessageBox.critical(self, "Error", f"Target '{dropdown_text}' not found in loaded data.")
                return

            file_path, self.current_target, self.current_redshift = self.target_data[matched_key]

        # Load the data for this target
            with open(file_path, 'rb') as f:
                self.spec_data = pickle.load(f)

            ignored_keys = {'Sightline', 'Target'}
            self.ion_list = [k for k in self.spec_data.keys() if k not in ignored_keys]

        # Process ion names
            self.ions = self.ion_name(self.ion_list)

            unique_ions = []
            [unique_ions.append(x) for x in self.ions if x not in unique_ions]

        # Populate ion dropdown
            self.ion_dropdown.clear()
            self.ion_dropdown.addItems(sorted(unique_ions))

        # Clear plot cache when target changes
            self.plot_cache = {}

        # Select first ion
            if unique_ions:
                self.ion_dropdown.setCurrentIndex(0)
                self.ion_changed()  # Simulate selection

            self.next_button.setEnabled(True)
            self.display_summary_table()

        except Exception as e:
            traceback.print_exc()
            qtw.QMessageBox.critical(self, "Error", f"Error loading data file: {str(e)}")

        finally:
        # Restore cursor
            qtw.QApplication.restoreOverrideCursor()
    
    def ion_name(self, list_of_ions):
        return ['NI' if ion.partition(' ')[0].startswith('NI') and
         len(ion.partition(' ')[0]) > 2 and
          len(ion.partition(' ')[0]) > 2 and (not(ion.partition(' ')[0].startswith('NII')))
            else ion.partition(' ')[0] for ion in list_of_ions]
    
    def ion_changed(self):
        """Handle ion selection change"""
        # Show busy cursor
        qtw.QApplication.setOverrideCursor(qtc.Qt.WaitCursor)
        
        try:
            selected = self.ion_dropdown.currentText()
            if not selected:
                return
                
            self.current_ion = selected
            
            # Update UI elements
            self.update_line_checkboxes()
            self.update_plot()
            self.update_button.setEnabled(True)
            
            # Check if this ion has already been measured
            self.check_existing_measurement()
            self.update_column_density()
        
        finally:
            # Restore cursor
            qtw.QApplication.restoreOverrideCursor()
        
    def checkbox_callback(self, trans, checkbox):
        self.spec_data[trans]['include_in_avg'] = checkbox.isChecked()
        self.line_selection_changed(trans)

    def find_valid_doublet_partners(self, trans, detected_ncols, max_delta=0.13, doublet_ratio=2.0, tol=0.1):
        """Find valid doublet partners for a given transition."""
        d1 = self.spec_data[trans]
        partners = []

        for tr2, nc2 in detected_ncols:
            if tr2 == trans:
                continue
            d2 = self.spec_data[tr2]
            try:
                f1 = float(d1['f'])
                f2 = float(d2['f'])
                wl1 = float(d1.get('lam_0', d1.get('wrest', 0)))
                wl2 = float(d2.get('lam_0', d2.get('wrest', 0)))
                r = max(f_lambda(wl1, f1), f_lambda(wl2, f2)) / min(f_lambda(wl1, f1), f_lambda(wl2, f2))
                delta_ncol = abs(d1['ncol_pyn'] - nc2)

                if abs(r - doublet_ratio) < tol and delta_ncol < max_delta:
                    partners.append((tr2, nc2, d2))
            except Exception:
                continue

        return partners

    def update_line_checkboxes(self):
        # Clear previous widgets
        for i in reversed(range(self.lines_scroll_layout.count())):
            child = self.lines_scroll_layout.itemAt(i).widget()
            if child:
                child.setParent(None)

        self.line_vars = {}
        self.quality_vars = {}
        self.entry_widgets = {}

        # Precompute detected ncols for current ion
        detected_ncols = []
        for j, trans in enumerate(self.ion_list):
            if self.ions[j] == self.current_ion:
                ion_data = self.spec_data[trans]
                if ion_data.get('nav_flag') in [0, -2] and np.isfinite(ion_data.get('ncol_pyn')):
                    detected_ncols.append((trans, ion_data['ncol_pyn']))

        for j, trans in enumerate(self.ion_list):
            if self.ions[j] != self.current_ion:
                continue

            ion_data = self.spec_data[trans]
            nav_flag = ion_data.get("nav_flag", np.nan)
            quality_flag = ion_data.get("quality_flag", "0")

            if nav_flag == -1:
                ncol = ion_data.get("ncol_linear2sig", np.nan)
                nerrl = np.nan
                nerrh = np.nan
            else:
                ncol = ion_data.get("ncol_pyn", np.nan)
                nerrl = ion_data.get("ncol_err_lo", np.nan)
                nerrh = ion_data.get("ncol_err_hi", np.nan)

            # Set row background color
            if nav_flag == 0:
                color = "#d8fdd8"
            elif nav_flag == -1:
                color = "#fffecb"
            else:
                color = "#fddede"

            # Row layout
            row_widget = qtw.QWidget()
            row_layout = qtw.QHBoxLayout(row_widget)
            row_widget.setFixedHeight(28)
            row_layout.setContentsMargins(1, 0, 1, 0)
            row_widget.setStyleSheet(f"background-color: {color};")

            # Transition label
            lbl_trans = qtw.QLabel(trans)
            lbl_trans.setMinimumWidth(120)
            lbl_trans.setAlignment(qtc.Qt.AlignLeft)
            row_layout.addWidget(lbl_trans)

            # Averaging checkbox
            default_checked = ion_data.get("include_in_avg", True)
            checkbox = qtw.QCheckBox()
            checkbox.setChecked(default_checked)
            self.spec_data[trans]['include_in_avg'] = checkbox.isChecked()
            checkbox.stateChanged.connect(lambda state, tr=trans, cb=checkbox: self.checkbox_callback(tr, cb))
            row_layout.addWidget(checkbox)
            self.line_vars[trans] = checkbox

            # ncol and error fields
            ncol_edit = qtw.QLineEdit(f"{ncol:.2f}" if np.isfinite(ncol) else "nan")
            ncol_edit.setReadOnly(True)
            ncol_edit.setMaximumWidth(80)
            ncol_edit.setAlignment(qtc.Qt.AlignCenter)
            row_layout.addWidget(ncol_edit)

            err_lo_edit = qtw.QLineEdit(f"{nerrl:.2f}" if np.isfinite(nerrl) else "nan")
            err_lo_edit.setReadOnly(True)
            err_lo_edit.setMaximumWidth(60)
            err_lo_edit.setAlignment(qtc.Qt.AlignCenter)
            row_layout.addWidget(err_lo_edit)

            err_hi_edit = qtw.QLineEdit(f"{nerrh:.2f}" if np.isfinite(nerrh) else "nan")
            err_hi_edit.setReadOnly(True)
            err_hi_edit.setMaximumWidth(60)
            err_hi_edit.setAlignment(qtc.Qt.AlignCenter)
            row_layout.addWidget(err_hi_edit)

            # nav_flag display
            flag_edit = qtw.QLineEdit(f"{int(nav_flag)}" if np.isfinite(nav_flag) else "nan")
            flag_edit.setReadOnly(True)
            flag_edit.setMaximumWidth(40)
            flag_edit.setAlignment(qtc.Qt.AlignCenter)
            row_layout.addWidget(flag_edit)

            # Quality dropdown
            # Analyst comment dropdown (separate from nav_flag or quality_flag)
            quality_combo = qtw.QComboBox()
            quality_combo.addItems(["0 - OK", "1 - Contaminated", "2 - Saturated"])

            # Load any previously saved comment index; fallback to 0
            try:
                current_comment = int(str(ion_data.get('comment', '0')).strip()[0])
            except:
                current_comment = 0

            quality_combo.setCurrentIndex(current_comment)
            quality_combo.setMaximumWidth(120)

            def on_quality_comment_changed(idx, tr=trans):
                self.spec_data[tr]['comment'] = str(idx)
                print(f"Updated comment for {tr} → {idx}")

            quality_combo.currentIndexChanged.connect(on_quality_comment_changed)

            row_layout.addWidget(quality_combo)
            self.quality_vars[trans] = quality_combo

            # Velocity display
            va = ion_data.get('va', np.nan)
            va_err = ion_data.get('va_err', np.nan)
            va_text = "—" if np.isnan(va) or np.isnan(va_err) else f"{va:.1f} ± {va_err:.1f} km/s"
            va_label = qtw.QLabel(va_text)
            va_label.setFixedWidth(110)
            va_label.setAlignment(qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
            row_layout.addWidget(va_label)

            # Show "Correct" button if valid doublet partner exists
            show_correct = False
            if nav_flag in [-2, 0] and np.isfinite(ncol) and ion_data.get("corrected_flag", 0) == 0:
                partners = self.find_valid_doublet_partners(trans, detected_ncols)
                if partners:
                    d2 = partners[0][2]
                    f1 = float(ion_data['f'])
                    f2 = float(d2['f'])

                    flag1 = nav_flag
                    flag2 = d2['nav_flag']


        # Ensure one is saturated, one is good
                    if (flag1 == -2 and flag2 == 0) or (flag1 == 0 and flag2 == -2):
            # Show the button ONLY on the weaker unsaturated line
                        if flag1 == 0 and f1 < f2:
                            show_correct = True
                            print("Show Correct button here (line is weak, flag 0)")
                        elif flag2 == 0 and f2 < f1:
                # Special case: you landed on the saturated line, but its partner is the weaker one
                # so do NOT show the button here
                            show_correct = False
                            print("Skip button: this is strong line, partner is weak")

            if show_correct:
                def make_correct_button(tr=trans, entry=ncol_edit):
                    def correct_ncol():
                        partners = self.find_valid_doublet_partners(tr, detected_ncols)
                        if not partners:
                            qtw.QMessageBox.warning(self, "Correction Skipped", "No valid doublet partner found.")
                            return

                        tr2, nc2, d2 = partners[0]

                        # Identify weak and strong
                        n1 = self.spec_data[tr]['ncol_pyn']
                        n2 = nc2
                        e1_lo = self.spec_data[tr].get('ncol_err_lo', 0)
                        e1_hi = self.spec_data[tr].get('ncol_err_hi', 0)
                        e2_lo = d2.get('ncol_err_lo', 0)
                        e2_hi = d2.get('ncol_err_hi', 0)

        # Compute symmetric errors
                        e1 = np.mean([e1_lo, e1_hi]) if np.isfinite(e1_lo) and np.isfinite(e1_hi) else 0
                        e2 = np.mean([e2_lo, e2_hi]) if np.isfinite(e2_lo) and np.isfinite(e2_hi) else 0

        # Determine weak vs strong
                        if n1 < n2:
                            logN = [n1, n2]
                            err_logN = [e1, e2]
                        else:
                            logN = [n2, n1]
                            err_logN = [e2, e1]

                        corrected_val, err_val = correct_saturation(logN, err_logN)

        # Show popup with correction proposal
                        original_val = self.spec_data[tr]['ncol_pyn']
                        original_err_lo = self.spec_data[tr].get('ncol_err_lo', 0)
                        original_err_hi = self.spec_data[tr].get('ncol_err_hi', 0)
        
                        msg = qtw.QMessageBox(self)
                        msg.setWindowTitle("Saturation Correction")
                        msg.setText(f"Correction calculated for {tr}:")
                        msg.setDetailedText(f"""
                Original value: {original_val:.2f} (+{original_err_hi:.2f}, -{original_err_lo:.2f})
                Corrected value: {corrected_val:.2f} ± {err_val:.2f}

                Partner transition: {tr2}
                Partner value: {nc2:.2f}
                """)
                        msg.setStandardButtons(qtw.QMessageBox.Yes | qtw.QMessageBox.No)
                        msg.setDefaultButton(qtw.QMessageBox.Yes)
        
                        reply = msg.exec_()
        
                        if reply == qtw.QMessageBox.Yes:
            # Store corrected values in separate fields - don't overwrite originals
                            self.spec_data[tr]['ncol_corrected'] = corrected_val
                            self.spec_data[tr]['ncol_corrected_err'] = err_val
                            self.spec_data[tr]['corrected_flag'] = 1
            
            # Add corrected label
                            corrected_lbl = qtw.QLabel(f"Corrected to {corrected_val:.2f}!")
                            corrected_lbl.setStyleSheet("color: green; font-weight: bold;")
                            row_layout.addWidget(corrected_lbl)
            
            # Update summary table and plot
                            self.update_column_density()
                            self.update_plot()
            
                            qtw.QMessageBox.information(self, "Success", 
                                f"Correction applied to {tr}. Press Update Measurement.")
                        else:
                            qtw.QMessageBox.information(self, "Cancelled", "Correction cancelled.")

                    return correct_ncol

                correct_btn = qtw.QPushButton("Correct")
                correct_btn.clicked.connect(make_correct_button())
                correct_btn.setMaximumWidth(80)
                row_layout.addWidget(correct_btn)

            elif ion_data.get('corrected_flag', 0) == 1:
                corrected_lbl = qtw.QLabel("Corrected!")
                corrected_lbl.setStyleSheet("color: green;")
                row_layout.addWidget(corrected_lbl)

            self.lines_scroll_layout.addWidget(row_widget)


    def update_plot(self):
        """Update the plot with current ion data based on selected lines"""
        self.plot_cache.clear()
        self.ax.clear()
    
        # Plot only transitions that are selected
        for j, trans in enumerate(self.ion_list):
            if self.ions[j] == self.current_ion and trans in self.line_vars and self.line_vars[trans].isChecked():
                # Check if we have cached data
                if trans not in self.plot_cache:
                    ion_data = self.spec_data[trans]
                
                # Safeguard against missing keys
                    if 'Nav' not in ion_data:
                        print(f"Warning: 'Nav' key missing for {trans}, skipping...")
                        continue
                
                    if 'vel' not in ion_data:
                        print(f"Warning: 'vel' key missing for {trans}, skipping...")
                        continue
                
                    if 'f' not in ion_data:
                        print(f"Warning: 'f' key missing for {trans}, skipping...")
                        continue
                
                    v1, v2 = -100, 100  # Default values if EWlims missing
                    if 'EWlims' in ion_data and ion_data['EWlims'] is not None:
                        v1, v2 = ion_data['EWlims']
                    
                    self.plot_cache[trans] = {
                        'vel': ion_data['vel'],
                        'Nav': ion_data['Nav'],
                        'f': ion_data['f'],
                        'v1': v1,
                        'v2': v2,
                        'contamination_mask': ion_data.get('contamination_mask', None)
                    }
            
            # Use cached data
                data = self.plot_cache[trans]
            
            # Additional safeguard: check if Nav data is valid
                if data['Nav'] is None or len(data['Nav']) == 0:
                    print(f"Warning: Invalid Nav data for {trans}, skipping...")
                    continue
            
            # Plot the main Nav profile
                self.ax.plot(data['vel'], data['Nav'], 
                        drawstyle='steps-mid', label=f"{trans}, f={data['f']}", zorder=3)
            
            # Replace the complex contamination shading section with this:
            # Add contamination shading if contamination_flag exists
                if data['contamination_mask'] is not None:
                    contamination_flag = np.array(data['contamination_mask'])
                    vel = np.array(data['vel'])
                    nav = np.array(data['Nav'])

                # Create a simple mask for contaminated regions
                    contaminated_mask = contamination_flag.astype(bool)

                    if np.any(contaminated_mask):
                        # Simple fill_between for all contaminated points
                        # Set non-contaminated points to NaN so they don't get filled
                        nav_masked = nav.copy()
                        nav_masked[~contaminated_mask] = np.nan
            
                        self.ax.fill_between(vel, nav_masked, alpha=0.3, color='red',
                                        step='mid', zorder=2)
                                    
            # Plot EW limit lines
                self.ax.axvline(data['v1'], linestyle='--', color='black', zorder=4)
                self.ax.axvline(data['v2'], linestyle='--', color='black', zorder=4)
    
    # Compute mean va and va_err from selected transitions
        va_vals = []
        va_err_vals = []
        for trans in self.line_vars:
            if self.line_vars[trans].isChecked():
                va = self.spec_data[trans].get('va')
                va_err = self.spec_data[trans].get('va_err')
                if va is not None and va_err is not None and np.isfinite(va) and np.isfinite(va_err):
                    va_vals.append(va)
                    va_err_vals.append(va_err)
    
        if va_vals:
            va_mean = np.mean(va_vals)
            va_err_mean = np.sqrt(np.mean(np.square(va_err_vals)))  # quadrature mean
        # Plot mean va as vertical line
            self.ax.axvline(va_mean, linestyle='--', color='#4400AA', zorder=2, alpha=0.5)
            self.ax.axvspan(va_mean - va_err_mean, va_mean + va_err_mean, 
                           color='#4400AA', alpha=0.07, zorder=1)
    
        self.ax.axhline(0, linestyle='--', color='k')
        self.ax.set_xlim(-500, 500)
        self.ax.set_xlabel('Velocity (km/s)')
        self.ax.set_ylabel('Na(v)')
        self.ax.set_title(f"{self.current_target}; z = {self.current_redshift}")
    
    # Only add legend if there are items to show
        if len(self.ax.get_lines()) > 1:  # More than just the horizontal line
            self.ax.legend()
    
        self.canvas.draw()

    def line_selection_changed(self, trans):
        """Handle line checkbox changes - update both plot and column density"""
    # Update column density calculation
        self.update_column_density()
    
    # Update the plot to reflect the new selection
        self.update_plot()

    def update_column_density(self):
        """Calculate and display column density, velocity, and redshift error based on selected lines"""
        N_values = []
        err_lo_values = []
        err_hi_values = []

        va_vals = []
        va_err_vals = []

        for trans, cb in self.line_vars.items():
            if cb.isChecked():
                spec = self.spec_data[trans].copy()
                # Use corrected values if they exist
                if spec.get('corrected_flag', 0) == 1:
                    spec['ncol_pyn'] = spec['ncol_corrected']
                    spec['ncol_err_lo'] = spec['ncol_corrected_err']
                    spec['ncol_err_hi'] = spec['ncol_corrected_err']
                elif trans in self.entry_widgets:
                # Use entry widget values if no stored correction exists
                    spec['ncol_pyn'] = self.entry_widgets[trans]['corrected_val']
                    spec['ncol_err_lo'] = self.entry_widgets[trans]['corrected_err']
                    spec['ncol_err_hi'] = self.entry_widgets[trans]['corrected_err']

                N_ret = self.column_density(spec)
                if np.isfinite(N_ret[0]):
                    N_values.append(N_ret[0])
                    err_lo_values.append(N_ret[1])
                    err_hi_values.append(N_ret[2])

            # Get velocity values
                va = spec.get('va')
                va_err = spec.get('va_err')
                if va is not None and va_err is not None and np.isfinite(va) and np.isfinite(va_err):
                    va_vals.append(va)
                    va_err_vals.append(va_err)

    # Compute average column density
        if N_values:
            if len(N_values) == 1:
                N_col = N_values[0]
                N_sig_low = err_lo_values[0]
                N_sig_high = err_hi_values[0]
            else:
                try:
                    if (np.any(np.isnan(err_hi_values)) or np.any(np.isnan(err_lo_values))):
                        linear_vals = 10 ** np.array(N_values)
                        mean_lin = np.mean(linear_vals)
                        N_col = np.log10(mean_lin)
                        N_sig_low = np.nan
                        N_sig_high = np.nan
                    else:
                        N_cal = logmean(np.array(N_values), np.array(err_hi_values), verbose=False)
                        N_col = N_cal[0]
                        N_sig_high = N_cal[1]
                        N_sig_low = -1.0 * N_cal[1]
                except Exception:
                    N_col, N_sig_low, N_sig_high = np.nan, np.nan, np.nan

            self.n_value_label.setText(f"{N_col:.3f}" if np.isfinite(N_col) else "nan")
            self.n_err_low_label.setText(f"{N_sig_low:.3f}" if np.isfinite(N_sig_low) else "nan")
            self.n_err_high_label.setText(f"{N_sig_high:.3f}" if np.isfinite(N_sig_high) else "nan")
        else:
            self.n_value_label.setText("nan")
            self.n_err_low_label.setText("nan")
            self.n_err_high_label.setText("nan")

    # Compute and display redshift error
        if va_err_vals:
            va_err_mean = np.sqrt(np.mean(np.square(va_err_vals)))  # quadrature mean
            c_kms = 299792.458
            z_err_mean = (va_err_mean / c_kms) * (1 + self.current_redshift)
            self.z_err_label.setText(f"{z_err_mean:.6f}")
        else:
            self.z_err_label.setText("nan")


    def column_density(self, spec):
        """Calculate column density from spectrum data"""
        N = 0
        err_N_lo = 0
        err_N_hi = 0
    
    # Follow the same calculation logic as the original code
        if ((spec['nav_flag']==0 or spec['nav_flag']==-2)):
            N = spec['ncol_pyn']
            err_N_lo = spec['ncol_err_lo']
            err_N_hi = spec['ncol_err_hi']
        elif (spec['nav_flag']==-1):
            N = spec['ncol_linear2sig']
            err_N_lo = np.nan
            err_N_hi = np.nan
        else:
            print("No Ncol!")
        
        return [N, err_N_lo, err_N_hi]

    def check_existing_measurement(self):
        """Check if this ion has already been measured for this target"""
        if self.summary is None:
            return

        for row in self.summary:
            if (row['targname'] == self.current_target and 
                row['ion'] == self.current_ion and 
                row['redshift'] == self.current_redshift):

                # Update UI with existing values
                self.n_value_label.setText(f"{row['N']:.3f}")
                self.n_err_low_label.setText(f"{row['N_sig_lo']:.3f}")
                self.n_err_high_label.setText(f"{row['N_sig_hi']:.3f}")

                # Set reliability and detection flags
                self.reliability_combo.setCurrentText(
                    f"{row['reliability']} - {self.get_reliability_text(row['reliability'])}"
                )
                self.detection_combo.setCurrentText(
                    f"{row['detection']} - {self.get_detection_text(row['detection'])}"
                )
                break


    def get_reliability_text(self, val):
        mapping = {
            1: "Multiple lines / Non-detection",
            2: "Atleast two lines / non-detection consistent",
            3: "Single line / multiple contaminations"
        }
        return mapping.get(val, "Unknown")

    def get_detection_text(self, val):
        mapping = {
            0: "Detection",
            -1: "Upper limit/Non-detection",
            -2: "Lower limit/Saturation"
        }
        return mapping.get(val, "Unknown")

    def update_measurement(self):
        """Update the measurement in the summary table quickly and reliably"""
        if self.current_ion == "":
            return

        try:
            N_col = float(self.n_value_label.text())
            N_sig_low = float(self.n_err_low_label.text())
            N_sig_high = float(self.n_err_high_label.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Invalid column density values")
            return

        selected_trans = [tr for tr, cb in self.line_vars.items() if cb.isChecked()]
        corrected_flag = 0

    # FIX: Check if ANY of the selected transitions have corrected values
        has_corrected = any(self.spec_data[tr].get('corrected_flag', 0) == 1 for tr in selected_trans)
        if has_corrected:
            corrected_flag = 1

    # Calculate velocity averages (this part is correct)
        va_vals = []
        va_err_vals = []
        for trans in selected_trans:
            va = self.spec_data[trans].get('va')
            va_err = self.spec_data[trans].get('va_err')
            if va is not None and va_err is not None and not np.isnan(va) and not np.isnan(va_err):
                va_vals.append(va)
                va_err_vals.append(va_err)

        if va_vals:
            va_mean = np.mean(va_vals)
            va_err_mean = np.sqrt(np.mean(np.square(va_err_vals)))
        else:
            va_mean = np.nan
            va_err_mean = np.nan
    
    # Calculate redshift error
        if np.isfinite(va_err_mean):
            c_kms = 299792.458
            z_err_mean = (va_err_mean / c_kms) * (1 + self.current_redshift)
        else:
            z_err_mean = np.nan

    # Get reliability and detection flags
        try:
            reliability = int(self.reliability_combo.currentText().split()[0])
            detection = int(self.detection_combo.currentText().split()[0])
        except Exception:
            reliability = 1
            detection = 0

    # FIX: Consistent quality flag handling
        quality_vals = []
        for tr in selected_trans:
            qbox = self.quality_vars.get(tr)
            if qbox:
                qval = qbox.currentIndex()  # Get index directly instead of parsing text
                quality_vals.append(str(qval))
                # Store consistently in spec_data
                self.spec_data[tr]['quality_flag'] = qval
                self.spec_data[tr]['comment'] = str(qval)  # Keep both for backward compatibility

        quality_flag = ";".join(quality_vals) if quality_vals else "0"

    # Update or add to summary table
        key = (self.current_target, self.current_ion, self.current_redshift)
        found = False
        for i, row in enumerate(self.summary):
            if (row['targname'], row['ion'], row['redshift']) == key:
                self.summary[i]['N'] = N_col
                self.summary[i]['N_sig_lo'] = N_sig_low
                self.summary[i]['N_sig_hi'] = N_sig_high
                self.summary[i]['va'] = va_mean
                self.summary[i]['va_err'] = va_err_mean
                self.summary[i]['z_err'] = z_err_mean
                self.summary[i]['reliability'] = reliability
                self.summary[i]['detection'] = detection
                self.summary[i]['corrected_flag'] = corrected_flag
                found = True
                break

        if not found:
            # Add new row
            new_row = [
                self.current_redshift, self.current_target, self.current_ion,
                va_mean, va_err_mean, z_err_mean, N_col, N_sig_low, N_sig_high,
                corrected_flag, reliability, detection
            ]
            self.summary.add_row(new_row)

    # Update the summary table display
        self.display_summary_table()
        QMessageBox.information(self, "Success", f"Measurement for {self.current_ion} updated")

    def next_ion(self):
        """Move to the next ion in the dropdown"""
        current_idx = self.ion_dropdown.currentIndex()
        total_items = self.ion_dropdown.count()

        if total_items == 0:
            QMessageBox.information(self, "Notice", "No ions to advance to.")
            return

        if current_idx < total_items - 1:
            self.ion_dropdown.setCurrentIndex(current_idx + 1)
            self.ion_changed()  # Trigger the change event
        else:
            QMessageBox.information(self, "Complete", "All ions processed for this target")

    def save_results(self):
        """Save the summary table and updated spec data to file"""
        if self.summary is None:
            return
    
        try:
            # Show busy cursor
            QApplication.setOverrideCursor(qtc.Qt.WaitCursor)

        # --- Save summary table ---
            self.summary.write('Summary/summary.csv', format='csv', overwrite=True)

        # --- Save current plot as PDF ---
            if self.current_ion and self.current_target:
                pdf_filename = f'Summary/{self.current_target}_z_{self.current_redshift}_summary_Nav.pdf'
                self.figure.savefig(pdf_filename)

            # --- Update quality_flag and corrected values ---
            for trans, qvar in self.quality_vars.items():
                self.spec_data[trans]['quality_flag'] = qvar.currentText()
    
        # Save updated pickle
        # First find the original path
            selected = self.target_dropdown.currentText()
            if selected in self.target_data:
                original_path, _, _ = self.target_data[selected]
                with open(original_path, 'wb') as f:
                    pickle.dump(self.spec_data, f)

            QMessageBox.information(self, "Success", "Results saved successfully")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving results: {str(e)}")
        finally:
        # Restore cursor
            QApplication.restoreOverrideCursor()

    def update_summary_table_entry(self, redshift, targ, ion, va, va_err, z_err, N_col, N_lo, N_hi, corr, rel, det):
        """Efficiently update one row of the summary table or insert it if new"""
        key = f"{targ}_{ion}_{redshift:.6f}"
        found = False

    # Search for existing row
        for row in range(self.summary_table.rowCount()):
            targ_item = self.summary_table.item(row, 1)
            ion_item = self.summary_table.item(row, 2)
            redshift_item = self.summary_table.item(row, 0)
        
            if targ_item and ion_item and redshift_item:
                row_key = f"{targ_item.text()}_{ion_item.text()}_{float(redshift_item.text()):.6f}"
                if row_key == key:
                # Update existing row
                    self.summary_table.item(row, 0).setText(f"{redshift:.6f}")
                    self.summary_table.item(row, 1).setText(targ)
                    self.summary_table.item(row, 2).setText(ion)
                    self.summary_table.item(row, 3).setText(f"{va:.3f}")
                    self.summary_table.item(row, 4).setText(f"{va_err:.3f}")
                    self.summary_table.item(row, 5).setText(f"{z_err:.2f}")
                    self.summary_table.item(row, 6).setText(f"{N_col:.3f}")
                    self.summary_table.item(row, 7).setText(f"{N_lo:.3f}")
                    self.summary_table.item(row, 8).setText(f"{N_hi:.3f}")
                    self.summary_table.item(row, 9).setText(str(corr))
                    self.summary_table.item(row, 10).setText(str(rel))
                    self.summary_table.item(row, 11).setText(str(det))
                    found = True
                    break

        if not found:
            # Add new row
            row_position = self.summary_table.rowCount()
            self.summary_table.insertRow(row_position)
        
            items = [
                f"{redshift:.6f}", targ, ion, f"{va:.3f}", f"{va:.3f}", f"{z_err:.2f}",
                f"{N_col:.3f}", f"{N_lo:.3f}", f"{N_hi:.3f}",
                str(corr), str(rel), str(det)
            ]
        
            for col, text in enumerate(items):
                item = QTableWidgetItem(text)
                self.summary_table.setItem(row_position, col, item)

    def display_summary_table(self):
        """Update the summary table widget to only show the current target/redshift"""
        if self.summary is None or self.current_target == "" or self.current_redshift == 0.0:
            return

        # Clear existing rows
        self.summary_table.setRowCount(0)

    # Add only rows matching the current target and redshift
        for row in self.summary:
            if row['targname'] != self.current_target or row['redshift'] != self.current_redshift:
                continue

            row_position = self.summary_table.rowCount()
            self.summary_table.insertRow(row_position)

            items = [
                f"{row['redshift']:.6f}",
                row['ion'],f"{row['va']:.1f}",
                f"{row['va_err']:.1f}",f"{row['z_err']:.2f}",
                f"{row['N']:.3f}",
                f"{row['N_sig_lo']:.3f}",
                f"{row['N_sig_hi']:.3f}",
                str(row['corrected_flag']),
                str(row['reliability']),
                str(row['detection']),
                    ]

            for col, text in enumerate(items):
                item = QTableWidgetItem(text)
                self.summary_table.setItem(row_position, col, item)
import math
class StackedPlotsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.setWindowTitle("Stacked Ion Profiles")
        self.setGeometry(100, 100, 1400, 900)
        self.setModal(False)
        # Create layout
        layout = QVBoxLayout()
        
        # Create tab widget
        self.tabs = qtw.QTabWidget()
        
        # Create tabs for different ion groups
        self.create_hi_tab()
        self.create_doublets_tab()
        self.create_all_tab()
        
        layout.addWidget(self.tabs)
        
        # Close button
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)
        
        self.setLayout(layout)
    
    def create_all_tab(self):
        """Create tab for all measured transitions (individual plots, not grouped)"""
        tab = QWidget()
        layout = QVBoxLayout()
    
    # Create scroll area for the plots
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
    
    # Get all individual transitions (not grouped)
        all_transitions = [trans for trans in self.parent.spec_data.keys() 
                          if trans not in ["Target", "Sightline"]]
    
        if all_transitions:
        # Create matplotlib figure with subplots
            n_plots = len(all_transitions)
            n_cols = 2
            n_rows = math.ceil(n_plots / n_cols)
        
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 2 * n_rows))
        
        # Store figure reference for saving
            self.all_transitions_fig = fig
        
        # Ensure axes is always a 2D array for consistent indexing
            if n_rows == 1 and n_cols == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            elif n_cols == 1:
                axes = axes.reshape(-1, 1)
        
        # Plot each transition individually
            for i, trans in enumerate(all_transitions):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col]
            
            # Plot as single transition (not as a group)
                self.plot_single_transition(ax, trans, plot_type='flux')
                
                # Only add x-axis label to bottom row
                if row == n_rows - 1:
                    ax.set_xlabel('Velocity (km/s)', fontsize=12)
                else:
                    ax.set_xlabel('')
                    ax.tick_params(labelbottom=False)
        
        # Hide unused subplots
            for i in range(n_plots, n_rows * n_cols):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col]
                ax.set_visible(False)
        
            plt.tight_layout(pad=0.5, h_pad=0.0, w_pad=0.3)
        
        # Add to scroll area
            canvas = FigureCanvas(fig)
            scroll_layout.addWidget(canvas)
        else:
        # No transitions found
            from PyQt5.QtWidgets import QLabel
            label = QLabel("No transitions found")
            scroll_layout.addWidget(label)
    
        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)
    
        tab.setLayout(layout)
        self.tabs.addTab(tab, "All Transitions")

    def create_hi_tab(self):
        """Create tab for HI transitions"""
        tab = QWidget()
        layout = QVBoxLayout()
    
    # Create scroll area for the plots
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
    
    # Get HI transitions and filter out non-spectral data
        hi_transitions_raw = self.get_transitions_for_ion("HI")
        hi_transitions = [trans for trans in hi_transitions_raw if trans not in ["Target", "Sightline"]]
    
        if hi_transitions:
        # Create matplotlib figure with subplots
            n_plots = len(hi_transitions)
            n_cols = 2
            n_rows = math.ceil(n_plots / n_cols)
        
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 2 * n_rows))
        
        # Store figure reference for saving
            self.hi_transitions_fig = fig
        
        # Ensure axes is always a 2D array for consistent indexing
            if n_rows == 1 and n_cols == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            elif n_cols == 1:
                axes = axes.reshape(-1, 1)
        
        # Plot each transition
            for i, trans in enumerate(hi_transitions):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col]
            
                self.plot_single_transition(ax, trans, plot_type='flux')
                
                # Only add x-axis label to bottom row
                if row == n_rows - 1:
                    ax.set_xlabel('Velocity (km/s)', fontsize=12)
                else:
                    ax.set_xlabel('')
                    ax.tick_params(labelbottom=False)
        
        # Hide unused subplots
            for i in range(n_plots, n_rows * n_cols):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col]
                ax.set_visible(False)
        
            plt.tight_layout(pad=0.5, h_pad=0.0, w_pad=0.3)
        
        # Add to scroll area
            canvas = FigureCanvas(fig)
            scroll_layout.addWidget(canvas)
        else:
            # No HI transitions found
            from PyQt5.QtWidgets import QLabel
            label = QLabel("No HI transitions found")
            scroll_layout.addWidget(label)
    
        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)
    
        tab.setLayout(layout)
        self.tabs.addTab(tab, "HI Transitions")

    def create_doublets_tab(self):
        """Create tab for doublet transitions (grouped in same plot)"""
        tab = QWidget()
        layout = QVBoxLayout()
    
    # Create scroll area for the plots
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
    
    # Get only the doublet groups (not singles)
        doublet_pairs = {
            'CIV': ['CIV 1548', 'CIV 1550'],
            'SiIV': ['SiIV 1393', 'SiIV 1402'],
            'OVI': ['OVI 1031', 'OVI 1037'],
            'MgII': ['MgII 2796', 'MgII 2803'],
            'NV': ['NV 1238', 'NV 1242'],
            'AlIII': ['AlIII 1854', 'AlIII 1862'],
            'OI+SiII': ['OI 1302', 'SiII 1304']  # User-defined doublet
        }
    
    # Get all available transitions
        all_transitions = [trans for trans in self.parent.spec_data.keys() 
                          if trans not in ["Target", "Sightline"]]
    
    # Get only doublet groups where at least one line is available
        doublet_groups = []
        for doublet_name, pair in doublet_pairs.items():
            available_pair = [trans for trans in pair if trans in all_transitions]
            if len(available_pair) >= 1:  # At least one line available
                doublet_groups.append(available_pair)
    
        if doublet_groups:
            # Create matplotlib figure with subplots
            n_plots = len(doublet_groups)
            n_cols = 2
            n_rows = math.ceil(n_plots / n_cols)
        
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 3 * n_rows))
            #fig.suptitle(f"Doublet Transitions - {self.parent.current_target}, z = {self.parent.current_redshift}", fontsize=14)
        
        # Store figure reference for saving
            self.doublets_fig = fig
        
        # Ensure axes is always a 2D array for consistent indexing
            if n_rows == 1 and n_cols == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            elif n_cols == 1:
                axes = axes.reshape(-1, 1)
        
        # Plot each doublet group
            for i, doublet_group in enumerate(doublet_groups):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col]
            
            # Plot the doublet group (both lines in same plot)
                self.plot_single_transition(ax, doublet_group, plot_type='Nav')
                
                # Only add x-axis label to bottom row
                if row == n_rows - 1:
                    ax.set_xlabel('Velocity (km/s)', fontsize=12)
                else:
                    ax.set_xlabel('')
                    ax.tick_params(labelbottom=False)
        
        # Hide unused subplots
            for i in range(n_plots, n_rows * n_cols):
                row = i // n_cols
                col = i % n_cols
                ax = axes[row, col]
                ax.set_visible(False)
        
            plt.tight_layout(pad=0.5, h_pad=0.0, w_pad=0.3)
        
        # Add to scroll area
            canvas = FigureCanvas(fig)
            scroll_layout.addWidget(canvas)
        else:
        # No doublet transitions found
            from PyQt5.QtWidgets import QLabel
            label = QLabel("No doublet transitions found")
            scroll_layout.addWidget(label)
    
        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)
    
        tab.setLayout(layout)
        self.tabs.addTab(tab, "Doublets")
    
    def get_transitions_for_ion(self, ion_name):
        """Get all transitions for a specific ion"""
        transitions = []
        for j, trans in enumerate(self.parent.ion_list):
            if self.parent.ions[j] == ion_name and trans in self.parent.spec_data:
                transitions.append(trans)
        return transitions

    def plot_single_transition(self, ax, trans, plot_type='Nav'):
        """Plot a single transition in the given axes
    
        Parameters:
        -----------
        ax : matplotlib axes object
        trans : str or list
            Single transition name or list of transitions (for doublets)
        plot_type : str
            'Nav' for Nav profile or 'flux' for normalized flux
        """
        # Check if ax is a valid matplotlib axes object
        if not hasattr(ax, 'plot'):
            print(f"Error: ax is not a matplotlib axes object. Type: {type(ax)}")
            return
    
    # Handle both single transitions and doublet lists
        if isinstance(trans, list):
            # Multiple transitions (doublets)
            transitions = trans
        else:
            # Single transition
            transitions = [trans]
    
        plotted_any = False
        legend_entries = []
    
        for trans_name in transitions:
            if trans_name not in self.parent.spec_data:
                continue
        
            ion_data = self.parent.spec_data[trans_name]
        
            # Add debug check to ensure ion_data is a dictionary
            if not isinstance(ion_data, dict):
                continue
        
        # Get plot data, handling potential missing keys
            vel = ion_data.get('vel')
            if plot_type == 'Nav':
                y_data = ion_data.get('Nav')
                y_label = 'Na(v)'
            else:  # plot_type == 'flux'
                y_data = ion_data.get('flux')
                y_label = 'Normalized Flux'
        
            f_val = ion_data.get('f')
        
            if vel is None or y_data is None or f_val is None:
                continue
        
        # Convert to numpy arrays for safety
            vel = np.array(vel)
            y_data = np.array(y_data)
        
        # Plot the main profile
            if len(transitions) == 1:
                # Single transition
                line, = ax.plot(vel, y_data, drawstyle='steps-mid', color='k', linewidth=1)
                legend_entries.append((line, f"{trans_name}, f={f_val:.3f}"))
            else:
                # Multiple transitions - use different colors
                colors = ['black', 'red', 'green', 'orange', 'purple', 'brown']
                color = colors[transitions.index(trans_name) % len(colors)]
                line, = ax.plot(vel, y_data, drawstyle='steps-mid', color=color, linewidth=1)
                legend_entries.append((line, f"{trans_name}, f={f_val:.3f}"))
        
        # Add contamination shading if available
            contamination_mask = ion_data.get('contamination_mask')
            if contamination_mask is not None:
                contamination_flag = np.array(contamination_mask)
                contaminated_mask = contamination_flag.astype(bool)
            
                if np.any(contaminated_mask):
                    y_masked = y_data.copy()
                    y_masked[~contaminated_mask] = np.nan
                    ax.fill_between(vel, y_masked, alpha=0.3, color='red', step='mid')
        
        # Add EW limits if available
            if 'EWlims' in ion_data and ion_data['EWlims'] is not None:
                v1, v2 = ion_data['EWlims']
                ax.axvline(v1, linestyle='--', color='b', alpha=0.7, linewidth=2)
                ax.axvline(v2, linestyle='--', color='b', alpha=0.7, linewidth=2)
        
        # Add va line if available
            va = ion_data.get('va')
            if va is not None and np.isfinite(va):
                ax.axvline(va, linestyle='--', color='purple', alpha=0.7, linewidth=2)
        
            plotted_any = True
    
        if not plotted_any:
            ax.text(0.5, 0.5, f"No data available", transform=ax.transAxes, 
                   ha='center', va='center')
    
    # Formatting
        if plot_type == 'Nav':
            ax.axhline(0, linestyle='--', color='g', alpha=0.5, linewidth=1)
        else:  # flux
            ax.axhline(1, linestyle='--', color='g', alpha=0.5, linewidth=1)
        if plot_type == 'flux':
            ax.axhline(0, linestyle='--', color='g', alpha=0.5, linewidth=1)
        
        ax.set_xlim(-500, 500)
        ax.set_ylabel(y_label, fontsize=12)
    
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3)
    
    # Add legend with ion names and f-values
        if legend_entries:
            lines, labels = zip(*legend_entries)
            ax.legend(lines, labels, fontsize=10, loc='upper right')

if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    
    # Create QApplication instance
    app = QApplication(sys.argv)
    
    # Create and show the main window
    window = SpectroscopicAnalysisGUI()
    window.show()
    
    # Start the event loop
    sys.exit(app.exec_())