// Include CsbScaffold
#load "../../.env/CsbScaffold.fsx"
#I @"../../.env/.aux/"
#load @"DeedleExtensions.fsx"
#nowarn "10001"
// If you want to use the wrappers for unmanaged LAPACK functions from of FSharp.Stats 
// include the path to the .lib folder manually to your PATH environment variable and make sure you set FSI to 64 bit

// use the following lines of code to ensure that LAPACK functionalities are enabled if you want to use them
// fails with "MKL service either not available, or not started" if lib folder is not included in PATH.
//open FSharp.Stats
//FSharp.Stats.Algebra.LinearAlgebra.Service()

open Deedle
open FSharpAux
open FSharp.Stats
open System.IO


let rawDataPath = Path.Combine [|__SOURCE_DIRECTORY__; "../" ;"data/PeptideRatioTable.tab";|]

/// Task 1: Read data & inspect the data Frame
let rawData :Frame<int,string> = 
    Frame.ReadCsv(rawDataPath,separators="\t") 

/// Task 2: Group peptides by ProtId
let groupedByProtID: Frame<string*int,string> = 
    rawData
    |> Frame.groupRowsBy "Protein"

/// Task 3: Get Columns with ratios
let ratioColumns = 
    groupedByProtID.ColumnKeys
    |> Seq.filter (fun x -> x.EndsWith("Log2Ratio"))

/// Task 5: Group peptides per protein name and compute the median ratio as a estimator for protein abundance.
//          This is done independently for technical replicates
let aggregatedPerProt (*:Frame<string*int,string>*) = 
    groupedByProtID
    |> Frame.sliceCols ratioColumns 
    |> Frame.getNumericCols
    |> Series.mapValues (Series.applyLevel Pair.get1Of2 (fun x -> x |> Series.values |> FSharp.Stats.Seq.median)) //*) |> Array.ofSeq |> FSharp.Stats.Quantile.OfSorted.compute 0.5))                   
    |> Frame.ofColumns

/// Task 6: According to the paper, exclude all proteins that have more than 10 missing values.
let filterProteinsWithMissingValues = 
    aggregatedPerProt
    |> Frame.sliceCols ratioColumns
    |> Frame.filterRowValues (fun x -> x.As<float>() |> Series.countValues >= 35)

/// task 7: Declare a function to calculate the mean row-wise. 
let calcMeansRowWise (columnKeys:string seq) (data: Frame<_,_>) = 
        data
        |> Frame.sliceCols columnKeys
        |> Frame.mapRowValues (fun row -> row.As<float>() |> Series.values |> FSharp.Stats.Seq.mean)

///TODO:
/// task 6: Declare a function to calculate the mean row-wise. 
let calcStdevRowWise (columnKeys:string seq) (data: Frame<_,_>) = 
        data
        |> Frame.sliceCols columnKeys
        |> Frame.mapRowValues (fun row -> row.As<float>() |> Series.values |> FSharp.Stats.Seq.stDev)

/// task 8: Create mapping from the biological replicate to the technical replicates.
let bioRepToTechnicalReplicates = 
    ratioColumns
    |> Seq.groupBy (fun x -> 
                        let tmp = x.Split('_')
                        tmp.[0]
                   ) 

/// task 9: Following the mapping from biological replicate to technical replicate,
///         Calculate the mean across technical replicates as an estimator for protein abundance
///         in a biological replicate.
let aggregateTechnicalReplicates = 
    bioRepToTechnicalReplicates 
    |> Seq.map (fun (bioRep, techRepColKeys) -> 
                   bioRep, calcMeansRowWise techRepColKeys filterProteinsWithMissingValues               
               )
    |> Frame.ofColumns

/// Task 10: Prepare a function to impute missing values.
let imp : (float [][] -> float [] [])= FSharp.Stats.ML.Impute.imputeBy (FSharp.Stats.ML.Impute.kNearestImpute 3) (Ops.isNan) //[|[|1.;1.1;1.2;nan|];[|1.;1.1;1.2;1.1|];[|1.;1.1;1.2;1.3|];[|1.;1.1;1.2;1.|];[|1.;1.1;1.2;nan|]|]

/// Task 11: Impute missing values. Leave Deedle for the first time
let withImputedValues =     
    let colKeys = aggregateTechnicalReplicates.ColumnKeys
    let rowKeys = aggregateTechnicalReplicates.RowKeys
    let dataForImputation = 
        aggregateTechnicalReplicates
        |> Frame.toArray2D 
        |> Array2D.toJaggedArray
    let imputedData = 
        imp dataForImputation    
        |> Array2D.ofJaggedArray
    Frame.ofArray2D imputedData
    |> Frame.indexColsWith colKeys
    |> Frame.indexRowsWith rowKeys

/// task 12: Create mapping from the timepoint to the biological replicates.
let timePointToBiologicalReplicates = 
    bioRepToTechnicalReplicates
    |> Seq.map fst
    |> Seq.groupBy (fun x -> 
                        let tmp = x.Split('-')
                        tmp.[0]
                   ) 


/// task 12: Following the mapping from biological replicate to technical replicate,
///         Calculate the mean across technical replicates as an estimator for protein abundance
///         in a biological replicate.
let meanAcrossBiologicalReplicates = 
    timePointToBiologicalReplicates 
    |> Seq.map (fun (bioRep, bioRepColKeys) -> 
                   bioRep + "_mean", calcMeansRowWise bioRepColKeys withImputedValues               
               )
    |> Frame.ofColumns

///
let stdevAcrossBiologicalReplicates =
    timePointToBiologicalReplicates 
    |> Seq.map (fun (bioRep, bioRepColKeys) -> 
                   bioRep + "_stdev", calcStdevRowWise bioRepColKeys withImputedValues               
               )
    |> Frame.ofColumns
///
let joinedF = Frame.mergeAll [withImputedValues;meanAcrossBiologicalReplicates;stdevAcrossBiologicalReplicates;] 

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
Chart.Column( FSharp.Stats.Distributions.Empirical.create bw testNAvsCA |> Map.toSeq )
|> Chart.Show

