# Compare drug sensitivities between two groups

S4 generic for comparing PRISM drug sensitivities between two groups
stored in a single \`OncoExperiment\` object.

## Usage

``` r
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

  An \`OncoExperiment\` object containing a top-level \`group\` column.

- group1:

  Label identifying group 1.

- group2:

  Label identifying group 2.

- unit:

  Character scalar describing the expected response unit.

- p_adj_method:

  Character scalar passed to \`stats::p.adjust()\`.

- effect_threshold:

  Non-negative numeric effect-size cutoff.

- ...:

  Additional method-specific arguments.

## Value

Method-specific comparison output.
