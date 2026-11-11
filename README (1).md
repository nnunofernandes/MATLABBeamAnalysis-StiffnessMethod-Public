# Beam Analysis in MATLAB Using the Stiffness Method

Teaching resources for solving continuous beams with the **direct stiffness method** in MATLAB. Two worked exercises connect the numerical implementation to analytical beam theory, including a case with support settlement.

The teaching package combines a PDF tutorial, guided MATLAB Live Script tutorials, starter scripts to complete during class, solved scripts, and detailed Word guides. Students can work through the calculation, inspect the intermediate matrices, and compare the numerical results with a hand solution using slope-deflection equations.

## Contents

- [Learning objectives](#learning-objectives)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Getting started](#getting-started)
- [Exercise 15–1 — Continuous beam under a point force and a distributed load](#exercise-151--continuous-beam-under-a-point-force-and-a-distributed-load)
- [Exercise 15–11 — Continuous beam with support settlement](#exercise-1511--continuous-beam-with-support-settlement)
- [Computational workflow](#computational-workflow)
- [Sign conventions and units](#sign-conventions-and-units)
- [Expected outputs](#expected-outputs)
- [Validation](#validation)
- [Suggested classroom activities](#suggested-classroom-activities)

## Learning objectives

After completing the exercises, students should be able to:

- Identify the degrees of freedom of an Euler–Bernoulli beam element.
- Construct element stiffness matrices and consistent nodal load vectors.
- Assemble a global stiffness matrix using element connectivity.
- Apply fixed supports, roller supports, and prescribed support displacements.
- Solve the reduced system for the unknown nodal rotation.
- Recover member end actions and support reactions.
- Calculate and interpret shear-force, bending-moment, and deflection diagrams.
- Compare the MATLAB solution with an analytical slope-deflection calculation.
- Verify equilibrium, compatibility, and linear superposition.

## Repository structure

The following layout defines the intended filenames for the complete teaching package. Place the PDF, starter scripts, and native `.mlx` versions under these names when adding them to the repository.

| Folder | File | Purpose |
|---|---|---|
| Root | `README.md` | Repository overview, exercise statements, and instructions. |
| `docs/` | `beam_stiffness_tutorial.pdf` | Main tutorial covering the theory and implementation. |
| `docs/` | `Exercise_15_1_MATLAB_and_Analytical_Guide.docx` | Detailed code explanation and analytical solution for Exercise 15–1. |
| `docs/` | `Exercise_15_11_MATLAB_and_Analytical_Guide.docx` | Detailed guide for Exercise 15–11, including settlement and superposition. |
| `tutorials/` | `exercise_15_1_beam_stiffness_live.mlx` | Guided Live Script tutorial for Exercise 15–1. |
| `tutorials/` | `exercise_15_11_beam_settlement_live.mlx` | Guided Live Script tutorial for Exercise 15–11. |
| `tutorials/` | `exercise_15_1_beam_stiffness_live_ready.m` | Commented source for creating the first Live Script in MATLAB. |
| `tutorials/` | `exercise_15_11_beam_settlement_live_ready.m` | Commented source for creating the settlement Live Script in MATLAB. |
| `matlab/starter/` | `exercise_15_1_starter.m` | **Before solution:** incomplete script for students to develop during class. |
| `matlab/starter/` | `exercise_15_11_starter.m` | **Before solution:** incomplete settlement exercise. |
| `matlab/solution/` | `exercise_15_1_beam_stiffness.m` | **After solution:** complete standalone solution for Exercise 15–1. |
| `matlab/solution/` | `exercise_15_11_beam_settlement.m` | **After solution:** complete standalone solution for Exercise 15–11. |

Starter files are intended for practice; solution files provide the completed calculation. The guided tutorials explain the purpose of each stage and the interpretation of its results.

## Requirements

- MATLAB for running the `.m` scripts.
- MATLAB Live Editor for the `.mlx` tutorials.
- A PDF viewer for the main tutorial.
- Microsoft Word or another compatible application for the optional `.docx` guides.

The completed MATLAB scripts use base MATLAB functionality and do not require additional toolboxes. They are standalone: each exercise can be run without the other exercise or an external data file.

## Getting started

1. Download or clone the repository and open its root folder in MATLAB.
2. Read the PDF tutorial in `docs/` to review the stiffness method and sign conventions.
3. Open the relevant starter script in `matlab/starter/` and attempt the exercise.
4. Work through the corresponding Live Script in `tutorials/`, running its sections in order.
5. Compare your implementation with the complete script in `matlab/solution/`.
6. Use the Word guide to follow the numerical calculation alongside the analytical solution.

To run a complete solution from the repository root:

```matlab
run(fullfile('matlab', 'solution', 'exercise_15_1_beam_stiffness.m'))
```

For the settlement exercise:

```matlab
run(fullfile('matlab', 'solution', 'exercise_15_11_beam_settlement.m'))
```

Run the exercises separately. Each script starts by clearing the workspace and closing existing figures. When using a tutorial, run its sections sequentially because later sections depend on variables created earlier.

### Creating a Live Script from its source

If you are using a `*_live_ready.m` tutorial source:

1. Open it in the MATLAB Editor.
2. Right-click its document tab and select **Open as Live Script**.
3. Save the new document as **MATLAB Live Code File (`*.mlx`)**, using the corresponding filename listed above.

MATLAB converts the publishing comments into formatted explanatory content. The source `.m` file can also be run directly.

## Exercise 15–1 — Continuous beam under a point force and a distributed load

### Problem statement

A continuous beam is fixed at supports **1** and **3** and rests on a roller at support **2**. The flexural rigidity **EI is constant** throughout the beam. Determine the moments at the two fixed supports.

| Item | Specification |
|---|---|
| Support 1 | Fixed, at global position `x = 0 m`. |
| Support 2 | Roller, at `x = 6 m`. |
| Support 3 | Fixed, at `x = 10 m`. |
| Span 1 | Length `6 m`; a `25 kN` downward point force acts at its midpoint, `x = 3 m`. |
| Span 2 | Length `4 m`; a downward uniform load of `15 kN/m` acts over the whole span. |
| Support movement | All prescribed support displacements are zero. |
| Flexural rigidity | Constant; its numerical value is unspecified. |

### Tasks

1. Represent the beam using three nodes and two beam elements.
2. Construct both element stiffness matrices.
3. Determine the consistent nodal loads for the midpoint force and the uniform load.
4. Assemble the global stiffness matrix and load vector.
5. Apply the supports and solve for the free rotation at node 2.
6. Recover the support reactions and member end moments.
7. Draw the shear-force, bending-moment, and stiffness-scaled deflection diagrams.
8. Compare the support moments with the analytical solution and check equilibrium.

**About EI:** with one common EI and zero imposed support movement, the reactions and internal moments are independent of its numerical value. The script uses a reference `EI = 1` for the force solution. Its raw rotations and deflections are reference values, not physical numerical predictions. The deflection plot shows **EI × v**; divide by the actual EI to obtain displacement in metres.

<details>
<summary>Show reference results for Exercise 15–1</summary>

Support reaction moments use counterclockwise positive:

| Support | Vertical reaction in kN | Reaction moment in kN·m | Internal bending moment in kN·m |
|---|---:|---:|---:|
| 1 | 12.375000 | +18.500000 | −18.500000 |
| 2 | 42.343750 | 0 | −19.250000 |
| 3 | 30.281250 | −20.375000 | −20.375000 |

The requested fixed-support reaction moments are **18.500 kN·m counterclockwise at support 1** and **20.375 kN·m clockwise at support 3**. Internal moments use sagging positive, so the support-section values are negative, indicating hogging.

</details>

## Exercise 15–11 — Continuous beam with support settlement

### Problem statement

A continuous beam is fixed at supports **1** and **3** and rests on a roller at support **2**. The intermediate support settles **30 mm downward**. The flexural rigidity is **EI = 3927 kN·m²**, constant throughout the beam. Determine the moments at the two fixed supports, including the effect of settlement.

| Item | Specification |
|---|---|
| Support 1 | Fixed, at global position `x = 0 m`. |
| Support 2 | Roller, at `x = 11 m`. |
| Support 3 | Fixed, at `x = 18 m`. |
| Span 1 | Length `11 m`; a downward uniform load of `30 kN/m` acts over the whole span. |
| Span 2 | Length `7 m`; a downward uniform load of `44 kN/m` acts over the whole span. |
| Support movement | Support 2 moves downward by `30 mm`, so `v₂ = −0.030 m`. |
| Flexural rigidity | `EI = 3927 kN·m²`. |

### Tasks

1. Form the element matrices and equivalent nodal load vectors.
2. Assemble the global stiffness matrix and load vector.
3. Prescribe the fixed-end conditions and the roller settlement.
4. Include the prescribed-displacement contribution in the reduced equation.
5. Solve for the rotation at node 2 and recover reactions and member end actions.
6. Plot shear force, bending moment, and physical deflection.
7. Separate the loading and settlement contributions using superposition.
8. Compare with the analytical slope-deflection solution and verify equilibrium and compatibility.

**About settlement:** impose the support movement as a known displacement, not as an additional force. Convert millimetres to metres before solving. Unlike the first exercise, the combined force solution depends on the specified EI because support movement produces stiffness-dependent internal actions.

<details>
<summary>Show reference results for Exercise 15–11</summary>

Support reaction moments use counterclockwise positive:

| Support | Vertical reaction in kN | Reaction moment in kN·m | Internal bending moment in kN·m |
|---|---:|---:|---:|
| 1 | 173.031244 | +333.895168 | −333.895168 |
| 2 | 324.056522 | 0 | −245.551481 |
| 3 | 140.912234 | −153.937116 | −153.937116 |

The unknown roller rotation is approximately **θ₂ = +0.035789006 rad**. The requested fixed-support reaction moments are **333.895 kN·m counterclockwise at support 1** and **153.937 kN·m clockwise at support 3**.

The reaction moments can also be checked by superposition:

| Support | Loads only in kN·m | Settlement only in kN·m | Combined in kN·m |
|---|---:|---:|---:|
| 1 | +326.384259 | +7.510909 | +333.895168 |
| 3 | −142.134259 | −11.802857 | −153.937116 |

</details>

## Computational workflow

Both exercises use the same global degree-of-freedom order:

$$
\mathbf{d}=\begin{bmatrix}v_1&\theta_1&v_2&\theta_2&v_3&\theta_3\end{bmatrix}^{T}
$$

Element 1 uses global indices `[1 2 3 4]`; element 2 uses `[3 4 5 6]`. The only free nodal degree of freedom is `θ₂`. The roller displacement is prescribed as zero in Exercise 15–1 and as `−0.030 m` in Exercise 15–11.

| Stage | Calculation | Main variables |
|---|---|---|
| Define the model | Enter geometry, connectivity, loading, EI, and settlement. | `xNode`, `connectivity`, `q`, `P`, `a`, `EI`, `settlement` |
| Number the unknowns | Map element degrees of freedom to global indices. | `elementDOF` |
| Form element equations | Calculate element stiffness and consistent load vectors. | `ke`, `fe` |
| Assemble | Add each element contribution at its global indices. | `K`, `F` |
| Apply constraints | Separate free and prescribed displacement entries. | `freeDOF`, `prescribedDOF`, `Kff`, `Kfc` |
| Solve | Determine the free rotation, including support movement. | `effectiveLoad`, `d` |
| Recover forces | Calculate reactions and member end actions. | `R`, `endAction`, `M_support` |
| Recover fields | Evaluate response throughout each span. | `xPlot`, `VPlot`, `MPlot`, `vPlot` |
| Verify | Check equilibrium, continuity, and reference results. | Residuals and `assert` statements |

The reduced stiffness equation is:

$$
\mathbf{K}_{ff}\mathbf{d}_f=\mathbf{F}_f-\mathbf{K}_{fc}\mathbf{d}_c
$$

The term involving the prescribed displacement vector is essential for the settlement exercise. MATLAB backslash solves the reduced system. Reactions and member end actions are then recovered from:

$$
\mathbf{R}=\mathbf{K}\mathbf{d}-\mathbf{F},\qquad
\mathbf{s}_e=\mathbf{k}_e\mathbf{d}_e-\mathbf{f}_e
$$

## Sign conventions and units

| Quantity | Positive convention | Unit |
|---|---|---|
| Vertical displacement and force | Upward | m and kN |
| Nodal rotation and reaction moment | Counterclockwise | rad and kN·m |
| Uniform transverse load | Upward; therefore the supplied downward loads are negative | kN/m |
| Internal bending moment | Sagging | kN·m |
| Shear force | Consistent with `V = dM/dx` in the diagrams | kN |
| Flexural rigidity EI | Positive stiffness | kN·m² |

A roller leaves rotation free and cannot supply an external reaction couple. **It does not create an internal hinge in a continuous beam.** The internal bending moment over that support can therefore be nonzero.

The member end-action convention and the internal diagram convention differ at the left end of an element: `M(0) = −Mi`, while `M(L) = Mj`. Check which convention is being used before comparing a reaction table with a bending-moment diagram.

## Expected outputs

Running a complete solution displays the global stiffness matrix and load vector, support reactions, internal support moments, and member end actions. It also produces a figure containing:

1. Beam geometry, support conditions, and applied loads.
2. Deflection or stiffness-scaled deflection.
3. Shear-force diagram.
4. Bending-moment diagram.

The settlement solution additionally displays the separate reaction-moment contributions from loading and support movement.

## Validation

The scripts include checks for:

- Symmetry of the global stiffness matrix.
- Whole-beam vertical-force and moment equilibrium.
- Zero reaction moment at the free roller rotation.
- Balance of the two member end moments at node 2.
- Agreement between recovered end displacements and the solved nodal values.
- Agreement with the reference support moments from slope-deflection calculations.
- Superposition of loading and settlement in Exercise 15–11.

When deliberately changing the problem data, disable the comparison with the original `referenceMoments` values. Keep the equilibrium and compatibility checks. Rounded values in this README may leave small residuals when substituted by hand.

## Suggested classroom activities

- Complete the starter files before viewing the solutions.
- Derive the equivalent nodal load vectors by hand and compare them with `F`.
- Explain why the beam has six global degree-of-freedom entries but only one unknown nodal displacement.
- Multiply the common EI in Exercise 15–1 by ten and predict which results change.
- Set the settlement in Exercise 15–11 to zero and recover the loading-only solution.
- Reverse the settlement sign and predict the change in the settlement-only actions.
- Explain why the shear diagram jumps at the roller while the bending-moment diagram remains continuous.

The exercises use linear Euler–Bernoulli beam theory. They do not include shear deformation, axial response, plasticity, or geometric nonlinearity.
