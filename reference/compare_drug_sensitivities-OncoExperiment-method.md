# Compare drug sensitivities between two OncoExperiment groups

Compares PRISM drug response values between two user-defined groups
stored in a single \`OncoExperiment\` object using per-drug Welch
two-sample t-tests.

## Usage

``` r
# S4 method for class 'OncoExperiment'
compare_drug_sensitivities(
  object,
  group1,
  group2,
  unit = NULL,
  p_adj_method = "BH",
  effect_threshold = 0,
  ...
)
```

## Arguments

- object:

  An \`OncoExperiment\` object with a top-level \`group\` column.

- group1:

  Label for the first group of models/samples.

- group2:

  Label for the second group of models/samples.

- unit:

  Optional character scalar describing the expected response unit (for
  example, \`"LFC"\`). If \`NULL\`, the unit is inferred from the assay
  metadata when available and normalized to \`"LFC"\` when the metadata
  uses an equivalent label such as \`"logfoldchange"\`.

- p_adj_method:

  Method passed to \`stats::p.adjust()\`.

- effect_threshold:

  Non-negative numeric threshold for \`abs(mean_diff)\` when computing
  significance flags.

- ...:

  Additional arguments. Currently unused.

## Value

A \`DrugSensitivityComparison\` object.

## Details

\*\*User is expected to define the groups prior to calling this
method.\*\* For example, to compare drug sensitivities between two
cancer types, subset the \`OncoExperiment\` and then assign groups:
\`COAD\$group \<- c(1, 1, 2, 1, 2, 2, 1, ...)\`

Whatever labels were chosen for the groups must be passed into the
function. The method looks in the top-level \`\$group\` column and
compares the selected values within the PRISM assay.

This method expects a loadable \`DrugResponseAssay\` (typically named
"PRISM") with assay data named "response".

The sign convention is: \`mean_diff = mean(group1) - mean(group2)\`.
Negative values indicate lower average response in \`group1\`.
