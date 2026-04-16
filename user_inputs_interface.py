import json
import os
import shutil
import sys
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Extraction scenarios folder creation
SCENARIO_FOLDER = "scenarios_save"

# Delete entire folder if it exists
if os.path.exists(SCENARIO_FOLDER):
    shutil.rmtree(SCENARIO_FOLDER)

# Recreate empty folder
os.makedirs(SCENARIO_FOLDER, exist_ok=True)

# Json file path (file will be created later when saving)
SCENARIO_FILE = os.path.join(SCENARIO_FOLDER, "user_inputs_config.json")

# JSON structure template (values will be lists in the file)
TEMPLATE = {
    "EXTRACTION": {
        "power_production_MW": [],
        "extraction_depth_m": [],
        "extraction_range_m": [],
        "extraction_water_m3s": [],
        "ch4_extraction_efficiency_percent": []
    },
    "REINJECTION": {
        "reinjection_depth_m": [],
        "reinjection_water_m3s": [],
        "ch4_reinjection_percent": [],
        "co2_reinjection_percent": []
    },
    "WASHING_EXTRACTION": {
        "washing_extraction_depth_m": [],
        "washing_extraction_range_m": [],
        "washing_extraction_water_m3s": []
    },
    "WASHING_REINJECTION": {
        "washing_reinjection_depth_m": [],
        "washing_reinjection_water_m3s": [],
        "ch4_washing_reinjection_percent": [],
        "co2_washing_reinjection_percent": []
    },
    "EXTRACTION_PERIOD": {
        "extraction_start_date": [],
        "extraction_end_date": []
    },
    "SIMULATION_MODEL": {
        "kivu_simstrat_path": []
    }
}

# Units for each parameter
UNITS = {
    "power_production_MW": "MW",
    "extraction_depth_m": "m",
    "extraction_range_m": "m",
    "extraction_water_m3s": "m³/s",
    "ch4_extraction_efficiency_percent": "%",

    "reinjection_depth_m": "m",
    "reinjection_water_m3s": "m³/s",
    "ch4_reinjection_percent": "%",
    "co2_reinjection_percent": "%",

    "washing_extraction_depth_m": "m",
    "washing_extraction_range_m": "m",
    "washing_extraction_water_m3s": "m³/s",

    "washing_reinjection_depth_m": "m",
    "washing_reinjection_water_m3s": "m³/s",
    "ch4_washing_reinjection_percent": "%",
    "co2_washing_reinjection_percent": "%",

    "extraction_start_date": "DD-MM-YYYY",
    "extraction_end_date": "DD-MM-YYYY",
    "kivu_simstrat_path": ""
}

class ScenarioUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Model Input Manager")

        self.data = self.load_or_init_data()
        self.scenario_count = self.get_existing_scenario_count()
        # Start on a new scenario (next index)
        self.current_scenario_index = self.scenario_count + 1

        self.fields = {}  # store widgets or widget tuples
        self.scenario_label_var = tk.StringVar()
        self.update_scenario_label()

        # Top-level layout: left column for scenario buttons, right for form
        top_container = tk.Frame(root)
        top_container.pack(side="top", fill="both", expand=True)

        # Left: scenario buttons column
        self.scenario_buttons_frame = tk.Frame(top_container, padx=10, pady=10)
        self.scenario_buttons_frame.pack(side="left", fill="y")

        tk.Label(self.scenario_buttons_frame, text="Scenarios",
                 font=("Arial", 10, "bold")).pack(anchor="nw")
        
        load_btn = tk.Button(
            self.scenario_buttons_frame,
            text="Load Inputs",
            command=self.load_inputs_from_file,
            width=12
        )
        load_btn.pack(anchor="nw", pady=5) # added for load existing json file

        self.scenario_buttons_inner = tk.Frame(self.scenario_buttons_frame)
        self.scenario_buttons_inner.pack(anchor="nw", pady=5)

        # Right: scenario label + form
        right_container = tk.Frame(top_container, padx=10, pady=10)
        right_container.pack(side="right", fill="both", expand=True)

        # Scenario label at top of right side
        label_frame = tk.Frame(right_container)
        label_frame.pack(side="top", fill="x")
        tk.Label(label_frame, textvariable=self.scenario_label_var,
                 font=("Arial", 14, "bold")).pack(anchor="center", pady=5)

        # Main input panel
        self.main_form_frame = tk.Frame(right_container)
        self.main_form_frame.pack(side="top", fill="both", expand=True)

        self.build_form(self.main_form_frame)

        # Bottom buttons (delete left, save right)
        bottom_frame = tk.Frame(root, padx=10, pady=10)
        bottom_frame.pack(side="bottom", fill="x")

        delete_btn = tk.Button(bottom_frame, text="Delete Scenario",
                               command=self.delete_scenario)
        delete_btn.pack(side="left", anchor="w")

        save_btn = tk.Button(bottom_frame, text="Save Scenario",
                             command=self.save_scenario)
        save_btn.pack(side="right", anchor="e")

        # Build scenario buttons for existing scenarios
        self.build_scenario_buttons()

    # ---------- Data handling ----------

    def pretty_label(self, key):
        # Remove units suffixes
        cleaned = key
        for suffix in ["_MW", "_m", "_m3s", "_percent"]:
            if cleaned.endswith(suffix):
                cleaned = cleaned.replace(suffix, "")
        # Replace underscores with spaces and capitalize
        cleaned = cleaned.replace("_", " ").strip().capitalize()
        return cleaned
    
    def load_inputs_from_file(self):
        # Ask user to select a JSON file
        filepath = filedialog.askopenfilename(
            title="Select Scenario JSON File",
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")]
        )
        if not filepath:
            return

        # Load the external JSON
        try:
            with open(filepath, "r") as f:
                imported = json.load(f)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{e}")
            return

        # Determine number of scenarios in the imported file
        # Use any parameter list length
        imported_count = None
        for section, items in imported.items():
            for key, lst in items.items():
                imported_count = len(lst)
                break
            if imported_count is not None:
                break

        if imported_count is None:
            messagebox.showerror("Error", "Invalid scenario file format.")
            return

        # Append imported scenarios to existing ones
        for section, items in TEMPLATE.items():
            for key in items:
                existing_list = self.data[section][key]
                imported_list = imported[section][key]

                # Append each imported scenario value
                for val in imported_list:
                    existing_list.append(val)

        # Update scenario count
        self.scenario_count += imported_count

        # Rebuild scenario buttons
        self.build_scenario_buttons()

        # Save updated data to config.json
        with open(SCENARIO_FILE, "w") as f:
            json.dump(self.data, f, indent=2)

        # Load the last imported scenario into UI
        self.current_scenario_index = self.scenario_count
        self.update_scenario_label()
        self.load_scenario(self.current_scenario_index)

        messagebox.showinfo("Success", f"Loaded {imported_count} scenarios.")


    def load_or_init_data(self):
        if os.path.exists(SCENARIO_FILE):
            with open(SCENARIO_FILE) as f:
                return json.load(f)
        # Initialize empty lists according to TEMPLATE
        data = {}
        for section, items in TEMPLATE.items():
            data[section] = {}
            for key in items:
                data[section][key] = []
        return data

    def get_existing_scenario_count(self):
        # Use any parameter list length as scenario count (if any exist)
        for section, items in self.data.items():
            for key, lst in items.items():
                return len(lst)
        return 0

    def update_scenario_label(self):
        self.scenario_label_var.set(f"Scenario {self.current_scenario_index}")

    # ---------- UI building ----------

    def build_form(self, parent):
        for section, items in TEMPLATE.items():
            # Skip SIMULATION_MODEL in main window
            if section == "SIMULATION_MODEL":
                continue

            frame = tk.LabelFrame(parent, text=section, padx=10, pady=10)
            frame.pack(fill="x", pady=5)

            for key in items:
                row = tk.Frame(frame)
                row.pack(fill="x", pady=1)

                #tk.Label(row, text=key, width=35, anchor="w").pack(side="left")
                tk.Label(row, text=self.pretty_label(key), width=35, anchor="w").pack(side="left")


                field_key = f"{section}.{key}"

                # Date fields: three numeric boxes (day, month, year)
                if key in ("extraction_start_date", "extraction_end_date"):
                    day_var = tk.StringVar(value="1")
                    month_var = tk.StringVar(value="1")
                    year_var = tk.StringVar(value="2016")

                    day_spin = tk.Spinbox(row, from_=1, to=31, width=3,
                                          textvariable=day_var)
                    month_spin = tk.Spinbox(row, from_=1, to=12, width=3,
                                            textvariable=month_var)
                    year_spin = tk.Spinbox(row, from_=1678, to=2200, width=5,
                                           textvariable=year_var)

                    day_spin.pack(side="left")
                    tk.Label(row, text="-").pack(side="left")
                    month_spin.pack(side="left")
                    tk.Label(row, text="-").pack(side="left")
                    year_spin.pack(side="left")

                    # Store tuple with a marker for date
                    self.fields[field_key] = ("date", day_spin, month_spin, year_spin)

                    # Delete button clears all three
                    #tk.Button(
                        #row,
                        #text="X",
                        #command=lambda ds=day_spin, ms=month_spin, ys=year_spin: self.clear_date(ds, ms, ys)
                    #).pack(side="left", padx=2)
                else:
                    entry = tk.Entry(row, width=20)
                    entry.pack(side="left")
                    self.fields[field_key] = ("entry", entry)
                    # --- deleting command on left pars
                    #tk.Button(
                        #row,
                        #text="X",
                        #command=lambda e=entry: e.delete(0, tk.END)
                    #).pack(side="left", padx=2)

                # Units label
                unit_text = UNITS.get(key, "")
                tk.Label(row, text=unit_text, width=12, anchor="w").pack(side="left")

    def build_scenario_buttons(self):
        # Clear existing buttons
        for child in self.scenario_buttons_inner.winfo_children():
            child.destroy()

        for i in range(1, self.scenario_count + 1):
            btn = tk.Button(
                self.scenario_buttons_inner,
                text=f"Scenario {i}",
                width=12,
                command=lambda idx=i: self.load_scenario(idx)
            )
            btn.pack(anchor="nw", pady=1)

    # ---------- Helpers ----------

    def clear_date(self, day_spin, month_spin, year_spin):
        day_spin.delete(0, tk.END)
        month_spin.delete(0, tk.END)
        year_spin.delete(0, tk.END)

    def clear_all_fields(self):
        for field_key, info in self.fields.items():
            kind = info[0]
            if kind == "entry":
                entry = info[1]
                entry.delete(0, tk.END)
            elif kind == "date":
                _, ds, ms, ys = info
                self.clear_date(ds, ms, ys)

    def parse_number(self, s):
        s = s.strip()
        if not s:
            return None
        try:
            i = int(s)
            return i
        except ValueError:
            try:
                f = float(s)
                # If it's an integer-like float, store as int
                if f.is_integer():
                    return int(f)
                return f
            except ValueError:
                return None

    # ---------- Scenario loading ----------

    def load_scenario(self, index):
        self.current_scenario_index = index
        self.update_scenario_label()

        # Fill fields from data lists
        for section, items in TEMPLATE.items():
            if section == "SIMULATION_MODEL":
                continue
            for key in items:
                field_key = f"{section}.{key}"
                kind, *widgets = self.fields[field_key]

                lst = self.data[section][key]
                value = ""
                if len(lst) >= index:
                    value = lst[index - 1]

                if kind == "entry":
                    entry = widgets[0]
                    entry.delete(0, tk.END)
                    if value is not None and value != "":
                        entry.insert(0, str(value))
                else:  # date
                    ds, ms, ys = widgets
                    ds.delete(0, tk.END)
                    ms.delete(0, tk.END)
                    ys.delete(0, tk.END)
                    if isinstance(value, str) and value:
                        # Expect "DD-MM-YYYY"
                        try:
                            d_str, m_str, y_str = value.split("-")
                            ds.insert(0, d_str)
                            ms.insert(0, m_str)
                            ys.insert(0, y_str)
                        except ValueError:
                            pass

    # ---------- Button actions ----------

    def delete_scenario(self):
        # Just clear current inputs
        self.clear_all_fields()

    def save_scenario(self):
        # Collect values for this scenario index
        index = self.current_scenario_index
        existing = self.scenario_count  # before potential extension

        for section, items in TEMPLATE.items():
            if section == "SIMULATION_MODEL":
                continue
            for key in items:
                field_key = f"{section}.{key}"
                kind, *widgets = self.fields[field_key]

                if kind == "entry":
                    entry = widgets[0]
                    raw = entry.get()
                    value = self.parse_number(raw)
                else:  # date
                    ds, ms, ys = widgets
                    d = ds.get().strip()
                    m = ms.get().strip()
                    y = ys.get().strip()
                    if d and m and y:
                        try:
                            d_i = int(d)
                            m_i = int(m)
                            y_i = int(y)
                            value = f"{d_i:02d}-{m_i:02d}-{y_i:04d}"
                        except ValueError:
                            value = ""
                    else:
                        value = ""

                lst = self.data[section][key]
                if index <= len(lst):
                    lst[index - 1] = value
                else:
                    lst.append(value)

        # Ensure SIMULATION_MODEL path list exists
        if "SIMULATION_MODEL" not in self.data:
            self.data["SIMULATION_MODEL"] = {"kivu_simstrat_path": []}
        elif "kivu_simstrat_path" not in self.data["SIMULATION_MODEL"]:
            self.data["SIMULATION_MODEL"]["kivu_simstrat_path"] = []

        # Update scenario count if we just added a new one
        if index > existing:
            self.scenario_count = index
            self.build_scenario_buttons()

        with open(SCENARIO_FILE, "w") as f:
            json.dump(self.data, f, indent=2)

        # After saving, open scenario popup
        self.open_scenario_popup()


    def validate_all_scenarios(self):
        """
        Returns:
            (True, None) if all scenarios are complete
            (False, scenario_index) if a scenario is incomplete
        """
        scenario_count = self.scenario_count

        for i in range(1, scenario_count + 1):
            for section, items in TEMPLATE.items():
                if section == "SIMULATION_MODEL":
                    continue

                for key in items:
                    value_list = self.data[section][key]

                    # If scenario index exceeds list length → missing
                    if i > len(value_list):
                        return False, i

                    value = value_list[i - 1]

                    # Empty string or None → missing
                    if value in ("", None):
                        return False, i

        return True, None

    # ---------- Scenario popup ----------

    def open_scenario_popup(self):
        popup = tk.Toplevel(self.root)
        popup.title(f"Scenario {self.current_scenario_index}")
        popup.geometry("300x120")
        popup.transient(self.root)
        popup.grab_set()
        popup.lift()

        tk.Label(popup, text=f"Scenario {self.current_scenario_index} saved.",
                 font=("Arial", 11)).pack(pady=10)

        btn_frame = tk.Frame(popup)
        btn_frame.pack(fill="x", pady=10, padx=10)

        add_btn = tk.Button(btn_frame, text="Add Scenario",
                            command=lambda: self.add_scenario_from_popup(popup))
        add_btn.pack(side="left", expand=True, fill="x", padx=5)

        sim_btn = tk.Button(btn_frame, text="Simulate",
                            command=lambda: self.simulate_from_popup(popup))
        sim_btn.pack(side="right", expand=True, fill="x", padx=5)

    def add_scenario_from_popup(self, popup):
        popup.destroy()
        self.current_scenario_index = self.scenario_count + 1
        self.update_scenario_label()
        self.clear_all_fields()

    #def simulate_from_popup(self, popup):
        #popup.destroy()
        #self.open_simulation_path_window()
    
    def simulate_from_popup(self, popup):
        popup.destroy()

        # Validate all scenarios before simulation
        ok, bad_index = self.validate_all_scenarios()

        if not ok:
            messagebox.showerror(
                "Missing Input",
                f"Scenario {bad_index} is incomplete.\n"
                "Please fill in all parameters before simulation."
            )
            # Jump user to the incomplete scenario
            self.load_scenario(bad_index)
            return

        # If everything is complete → proceed
        self.open_simulation_path_window()


    # ---------- Simulation path window ----------

    def open_simulation_path_window(self):
        win = tk.Toplevel(self.root)
        win.title("Simulation Path")
        win.geometry("500x160")
        win.transient(self.root)
        win.grab_set()
        win.lift()

        tk.Label(win, text="Select Kivu Simstrat model path:",
                 font=("Arial", 11)).pack(pady=10)

        path_frame = tk.Frame(win)
        path_frame.pack(fill="x", padx=10)

        path_var = tk.StringVar()
        path_entry = tk.Entry(path_frame, textvariable=path_var, width=45)
        path_entry.pack(side="left", fill="x", expand=True)

        def browse():
            # Make sure dialog is on top of this window
            win.lift()
            path = filedialog.askdirectory(parent=win)
            if path:
                path_var.set(path)

        tk.Button(path_frame, text="Browse", command=browse).pack(side="left", padx=5)

        def start_simulation():
            path = path_var.get().strip()
            if not path:
                messagebox.showerror("Error", "Please select a model path first.")
                return

            run_script = os.path.join(path, "simulation_with_extraction_launcher.py")
            if not os.path.isfile(run_script):
                messagebox.showerror("Error", f"'simulation_with_extraction_launcher.py' not found in:\n{path}")
                return

            # Store path per scenario
            path_list = self.data["SIMULATION_MODEL"]["kivu_simstrat_path"]
            idx = self.current_scenario_index
            if idx <= len(path_list):
                path_list[idx - 1] = path
            else:
                # Fill missing with None if needed
                while len(path_list) < idx - 1:
                    path_list.append(None)
                path_list.append(path)

            with open(SCENARIO_FILE, "w") as f:
                json.dump(self.data, f, indent=2)

            # Run the model
            try:
                subprocess.Popen([sys.executable, run_script], cwd=path)
                messagebox.showinfo("Simulation", "Simulation started (simulation_with_extraction_launcher.py).")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to start simulation:\n{e}")

            win.destroy()

        tk.Button(win, text="Start Simulation", command=start_simulation).pack(pady=15)


if __name__ == "__main__":
    root = tk.Tk()
    app = ScenarioUI(root)
    root.mainloop()