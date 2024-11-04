from electrolyzer import Stack, Supervisor, run_electrolyzer_zbt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm # matplotlib's color map library
import matplotlib as mpl
import electrolyzer.inputs.validation as val
import plotly.graph_objects as go; import numpy as np
from plotly_resampler import FigureResampler, FigureWidgetResampler
import pandas as pd
import math

if __name__ == "__main__":
    # Initialize system
    fname_input_modeling = "./DEl_PM2_modeling_options.yaml"
    modeling_options = val.load_modeling_yaml(fname_input_modeling)
    elec_sys = Supervisor.from_dict(modeling_options["electrolyzer"])
    elec = elec_sys.stacks[0]
