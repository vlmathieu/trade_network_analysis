# Trade-profiles figure under alternative comparison windows (supplementary
# material). The main-text figure compares 2000 and 2020; these three show the
# same figure for the last pre-COVID year (2019), the first full post-shock
# year (2021), and the full study period (1996-2023), to show that the choice
# of comparison years does not drive the result.
#
# Bubble scale (size_max / size_breaks) is pinned to a common maximum across
# the three so they are comparable with each other. It differs from the
# main-text figure, whose own scale is left untouched — say so in the captions.
#
# Zoom-box overrides are per-window: import degrees shift enough between years
# to move countries in and out of the two inset boxes. 2000-2021 needs none
# (its qualifying set and box membership are identical to 2000-2020).
#
# Label nudges work the same way. nudge_main / nudge_a / nudge_b are merged on
# top of NUDGE_MAIN / NUDGE_A / NUDGE_B in the script, so a value set here moves
# a label in this figure alone and leaves the main-text and deflated figures
# untouched. Each entry is {"Country": [x, y]} in sqrt-transformed axis units.
# The inset badges are drawn at the top-left corner of their box and are opaque
# to ggrepel, so a label landing under one has to be nudged clear by hand.

PROFILE_WINDOW_SIZE_MAX    = 12000
PROFILE_WINDOW_SIZE_BREAKS = [1000, 6000, 12000]


rule plot_contributor_profiles_2000_2019:
    input:
        'results/network_analysis/agg_eu/output/contributor_profiles.csv'
    output:
        expand('results/network_analysis/agg_eu/plot/01/contributor_profiles_2000_2019.{ext}',
               ext = config['figure_ext'])
    params:
        fao_divisions = ['01'],
        year_start    = 2000,
        year_end      = 2019,
        ext           = config['figure_ext'],
        size          = 'primary_value',
        threshold     = config['threshold_main_contributors'],
        out_stem      = 'contributor_profiles_2000_2019',
        size_max      = PROFILE_WINDOW_SIZE_MAX,
        size_breaks   = PROFILE_WINDOW_SIZE_BREAKS,
        # New Zealand (14), Australia (14), Russia (13) and Canada (10) all sit
        # above the default ceiling of 10 import partners in 2019 and would
        # fall out of inset B into the crowded middle of the main panel.
        # 15 keeps the same eight-country inset as 2000-2020 and still leaves
        # Malaysia (16) out.
        zoom_a_y      = [25, 90],
        zoom_b_x      = [13, 36],
        zoom_b_y      = [0.5, 23.5],
        # The B badge sits at (10, 15) here and the Malaysia label runs under
        # it: push Malaysia clear.
        nudge_main    = {},
        nudge_a       = {"Rep. of Korea":   [ 0.70,  0.40],
                         "Hong Kong":       [-0.80,  0],
                         "Türkiye":         [-0.30,  0.30]},
        # Australia (22, 14), New Zealand (24, 14), Russia (27, 13) and Canada
        # (26, 10) bunch along the top edge of the inset.
        nudge_b       = {"Switzerland":     [-0.40,  0.10],
                         "Malaysia":        [ 0.00,  0.55],
                         "New Zealand":     [-0.25,  0.35],
                         "Russia":          [ 0.20,  0.25],
                         "Canada":          [ 0.25, -0.25],
                         "Australia":       [-0.20, -0.25],
                         "Cameroon":        [-0.20,  0.25]}
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_contributor_profiles.R'


rule plot_contributor_profiles_2000_2021:
    input:
        'results/network_analysis/agg_eu/output/contributor_profiles.csv'
    output:
        expand('results/network_analysis/agg_eu/plot/01/contributor_profiles_2000_2021.{ext}',
               ext = config['figure_ext'])
    params:
        fao_divisions = ['01'],
        year_start    = 2000,
        year_end      = 2021,
        ext           = config['figure_ext'],
        size          = 'primary_value',
        threshold     = config['threshold_main_contributors'],
        out_stem      = 'contributor_profiles_2000_2021',
        size_max      = PROFILE_WINDOW_SIZE_MAX,
        size_breaks   = PROFILE_WINDOW_SIZE_BREAKS,
        # Malaysia and Canada cross
        nudge_main    = {"Canada":      [-0.50, -0.50],
                         "Malaysia":    [ 0.50, -0.60]},
        # "New Zealand" (19, 9) and "Australia" (25, 9) overprint each other.
        nudge_b       = {"New Zealand": [-0.25, -0.25],
                         "Australia":   [-0.15, -0.15],
                         "Russia":      [-0.20, -0.20],
                         "Brazil":      [-0.25,  0.15]}
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_contributor_profiles.R'


rule plot_contributor_profiles_1996_2023:
    input:
        'results/network_analysis/agg_eu/output/contributor_profiles.csv'
    output:
        expand('results/network_analysis/agg_eu/plot/01/contributor_profiles_1996_2023.{ext}',
               ext = config['figure_ext'])
    params:
        fao_divisions = ['01'],
        year_start    = 1996,
        year_end      = 2023,
        ext           = config['figure_ext'],
        size          = 'primary_value',
        threshold     = config['threshold_main_contributors'],
        out_stem      = 'contributor_profiles_1996_2023',
        size_max      = PROFILE_WINDOW_SIZE_MAX,
        size_breaks   = PROFILE_WINDOW_SIZE_BREAKS,
        # Nothing reaches 80 import partners in either year, so the default
        # ceiling of 90 leaves the top third of inset A empty. The lower bound
        # stays at 27: dropping it would pull Canada (20 import partners in
        # 1996) into an inset captioned as the main Asian importers.
        zoom_a_y      = [26, 82],
        # New Zealand sits exactly on 10 import partners in 1996, so raise the
        # ceiling slightly. 11 keeps Malaysia (12 in 2023) in the main panel,
        # where its 1996 position also sits.
        zoom_b_x      = [9, 36],
        zoom_b_y      = [0.5, 11.5],
        # Three fixes, all in the main panel: the Malaysia label runs under the
        # B badge at (10, 11); Philippines, Thailand and Türkiye stack below
        # box A; Canada overlaps Switzerland. Inset B renders clean.
        nudge_main    = {"Türkiye":     [ 0.80, -0.05],
                         "Thailand":    [-0.90,  0],
                         "Philippines": [-0.50, -0.50],
                         "Norway":      [-0.90,  0.25],
                         "Canada":      [-0.90,  0.30],
                         "Malaysia":    [ 0.60,  0.40],
                         "Switzerland": [-1.20,  0]},
        nudge_a       = {"Rep. of Korea":   [-0.30,  0.30],
                         "Viet Nam":        [ 0.20,  0.80]}
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script:
        '../scripts/plot_contributor_profiles.R'
