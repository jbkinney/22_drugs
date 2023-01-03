```
README.md       # This file

data/
    molecular_dynamics_simulations/
        22.12.18_branaplam_enol.pdb   # Fig. S13 B (PDB file)
        22.12.18_branaplam_enol.xtc   # Fig. S13 B (MD trajectory)
        22.12.18_risdiplam.pdb        # Fig. S13 A (PDB file)
        22.12.18_risdiplam.xtc        # Fig. S13 A (MD trajectory)
    mpsa/
        psi_smn2_dmso.csv       # Fig. 1, 3, S1, S2, S3, S5, S11, S12
        psi_smn2_rg.csv         # Fig. 1, 3, S1, S2, S3, S5, S11, S12
        psi_smn2_nvs.csv        # Fig. 1, 3, S1, S2, S3, S5, S11, S12
        psi_elp1_dmso.csv       # Fig. S4
        psi_elp1_rg.csv         # Fig. S4
        psi_elp1_nvs.csv        # Fig. S4
    qPCR/
        dose_response_curves.xlsx       # Fig. 5E-L, 6A-K
        linear_mixture_curves.xlsx      # Fig. 6L-O
        smn1_345mutants.xlsx            # Fig. 4F,G
    rna-seq/
        allelic_manifold_fits.csv.gz    # Fig. 2D,E (fit curves)
        rmats_results.csv.gz            # Fig. 2E (data points)

fig1/
    fig1D-I_S1_S2_S3_S5_S8.ipynb        # SMN2 MPSA analysis
    figS4.ipynb                         # ELP1 MPSA analysis

fig2/
    fig2BDE_S7.ipynb   # allelic manifolds, effect scatter, hyp min* motif

fig3/
    fig3_S11_S12_inference.ipynb    # Do joint model inference
    fig3C-H_S12.ipynb                  # Plot Fig. 3C-H and S12
    fig3B_S11.ipynb                    # Plot Fig. 3B and S11

fig4/
    drugs.icl       # Fig. 4C,D,E (planar structures)
    figS13.ipynb    # dihedral histograms
    fig4FG.ipynb    # qPCR validation
    fig4A.pse       # U1/5'ss apo, PyMOL scene
    fig4B.pse       # U1/5'ss/SMN-C5 complex, PyMOL scene

fig5/
    fig5E-L_6A-K_inference.ipynb  # Dose-response curve inference, Figs. 5, 6 
    fig5D-L.ipynb                 # Fig. 5 dose-response curves
    all_dose_response_curves/     # Diagnostic plots of dose-response curves

fig6/
    fig6L-O_inference.ipynb       # Linear-mixture curve inference, Fig. 6
    fig6.ipynb                    # Fig. 6 dose-response & linear-mixture curves
    all_lienar_mixture_curves/    # Diagnostic plots of linear-mixutre curves

mcmc_samples/
    mcmc_*.pkl      # Bayesian posterior samples computed using numpyro.
```
