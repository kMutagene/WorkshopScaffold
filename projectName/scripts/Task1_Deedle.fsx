// Include CsbScaffold
#load "../../.env/CsbScaffold.fsx"
#I @"../../.env/.aux/"
// If you want to use the wrappers for unmanaged LAPACK functions from of FSharp.Stats 
// include the path to the .lib folder manually to your PATH environment variable and make sure you set FSI to 64 bit

// use the following lines of code to ensure that LAPACK functionalities are enabled if you want to use them
// fails with "MKL service either not available, or not started" if lib folder is not included in PATH.
//open FSharp.Stats
//FSharp.Stats.Algebra.LinearAlgebra.Service()


open BioFSharp
open Deedle
open FSharpAux
open Deedle
open FSharp.Stats
open BioFSharp.Elements
open System.IO
#load @"DeedleExtensions.fsx"

let rawDataPath = Path.Combine [|__SOURCE_DIRECTORY__; "../" ;"data/PeptideRatioTable.tab";|]

///// Task 1: read data & inspect the data Frame
let rawData :Frame<int,string> = 
    Frame.ReadCsv(rawDataPath,separators="\t") 

///// Task 2: filter out peptides with quantScore < 0.7 
let filteredByQuantScore = 
    rawData
    // TODO: 

//// Task 3: Group peptides by ProtId
let groupedByProtID: Frame<string*int,string> = 
    rawData
    |> Frame.groupRowsBy "Protein"

/// Task 4: get Columns with ratios
let ratioColumns = 
    groupedByProtID.ColumnKeys
    |> Seq.filter (fun x -> x.EndsWith("Log2Ratio"))

/// Task 5: group peptides per protein name and compute the median ratio as a estimator for protein abundance.
//          This is done independently for technical replicates
let aggregatedPerProt (*:Frame<string*int,string>*) = 
    groupedByProtID
    |> Frame.sliceCols ratioColumns 
    |> Frame.getNumericCols
    |> Series.mapValues (Series.applyLevel Pair.get1Of2 (fun x -> x |> Series.values |> Seq.filter (fun x -> isNan x |> not) |> Seq.sort |> FSharp.Stats.Seq.median))                   
    |> Frame.ofColumns


/// task 6: Declare a function to calculate the mean row-wise. 
let calcMeans (columnKeys:string seq) (data: Frame<_,_>) = 
        data
        |> Frame.sliceCols columnKeys
        |> Frame.mapRowValues (fun row -> row.As<float>() |> Series.values |> FSharp.Stats.Seq.mean)


/// task 7: Create mapping from the biological replicate to the technical replicates.
let groupTechnicalReplicates = 
    ratioColumns
    |> Seq.groupBy (fun x -> 
                        let tmp = x.Split('_')
                        tmp.[0]
                   ) 


/// task 8: Following the mapping from biological replicate to technical replicate,
///         Calculate the mean across technical replicates as an estimator for protein abundance
///         in a biological replicate.
let aggregateTechnicalReplicates = 
    groupTechnicalReplicates 
    |> Seq.map (fun (bioRep, techRepColKeys) -> 
                   bioRep, calcMeans techRepColKeys aggregatedPerProt               
               )
    |> Frame.ofColumns


/// Task 9: Prepare a function to impute missing values.
let imp : (float [][] -> float [] [])= FSharp.Stats.ML.Impute.imputeBy (FSharp.Stats.ML.Impute.kNearestImpute 3) (Ops.isNan) //[|[|1.;1.1;1.2;nan|];[|1.;1.1;1.2;1.1|];[|1.;1.1;1.2;1.3|];[|1.;1.1;1.2;1.|];[|1.;1.1;1.2;nan|]|]

/// Task 10: Impute missing values. Leave Deedle for the first time
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
    |> Frame.mapColKeys (fun x -> colKeys |> Seq.item x)
    |> Frame.mapRowKeys (fun x -> rowKeys |> Seq.item x)


/// task 11: Create mapping from the timepoint to the biological replicates.
let groupBiologicalReplicates = 
    groupTechnicalReplicates
    |> Seq.map fst
    |> Seq.groupBy (fun x -> 
                        let tmp = x.Split('-')
                        tmp.[0]
                   ) 

/// task 6: Declare a function to calculate the mean row-wise. 
let calcstDev (columnKeys:string seq) (data: Frame<_,_>) = 
        data
        |> Frame.sliceCols columnKeys
        |> Frame.mapRowValues (fun row -> row.As<float>() |> Series.values |> FSharp.Stats.Seq.stDev)

/// task 12: Following the mapping from biological replicate to technical replicate,
///         Calculate the mean across technical replicates as an estimator for protein abundance
///         in a biological replicate.
let meanAcrossBiologicalReplicates = 
    groupBiologicalReplicates 
    |> Seq.map (fun (bioRep, bioRepColKeys) -> 
                   bioRep + "_mean", calcMeans bioRepColKeys withImputedValues               
               )
    |> Frame.ofColumns

///
let stdevAcrossBiologicalReplicates =
    groupBiologicalReplicates 
    |> Seq.map (fun (bioRep, bioRepColKeys) -> 
                   bioRep + "_stdev", calcstDev bioRepColKeys withImputedValues               
               )
    |> Frame.ofColumns

///
let joinedF = Frame.mergeAll [withImputedValues;meanAcrossBiologicalReplicates;stdevAcrossBiologicalReplicates;] 


//withImputedValues.SaveCsv(@"C:\Users\david\OneDrive - tukl\Dokumente\BioInfo\91_WorkshopBerlin\Interactive Deedling\data\f�rVenny.tab",separator='\t',includeRowKeys=true)                   
  
