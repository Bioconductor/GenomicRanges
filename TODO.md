# SummarizedExperiment

Long term
- Separate package
- Separate classes for 'DataFrame' and '\*Ranges' rowData

Immediate

1. Implement rowRanges() as a drop-in replacement for rowData()
2. Deprecate rowData() in favour of rowRanges()

Short term (start of next release cycle)

3. Implement separate classes for DataFrame and \*Ranges rowData()

Possibilities?

- SummarizedExperiment virtual base class with derived classes

  - SummarizedExperimentDF
    @rowData: DataFrame

  - SummarizedExperimentGR
    @rowRanges: *Ranges; rowData() == mcols(rowRanges())

- SummarizedExperiment as 'DataFrame' base class,
  SummarizedExperimentGR as derived class

  - SummarizedExperiment
    @rowData: DataFrame

  - SummarizedExperimentGR
    @rowRanges: \*Ranges; no mcols() (rowData() from inheritted
    @rowData)

