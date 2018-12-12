// Include CsbScaffold
#load "../../.env/CsbScaffold.fsx"
#load "Task_3_2_Deedle.fsx"
#nowarn "10001"
// If you want to use the wrappers for unmanaged LAPACK functions from of FSharp.Stats 
// include the path to the .lib folder manually to your PATH environment variable and make sure you set FSI to 64 bit

// use the following lines of code to ensure that LAPACK functionalities are enabled if you want to use them
// fails with "MKL service either not available, or not started" if lib folder is not included in PATH.
//open FSharp.Stats
//FSharp.Stats.Algebra.LinearAlgebra.Service()

open System
open FSharpAux
open Task_3_2_Deedle
open Deedle
open FSharp.Stats
open FSharpGephiStreamer
#time

let applyThreshold thr m = Array2D.map (fun x -> if x > thr then x else 0.) m

let timeSeriesData : Frame<string,string> = 
    let path = __SOURCE_DIRECTORY__ + @"..\..\data\AggregatedProteinQuantTable.tab"
    Frame.ReadCsv(path,indexCol = "Key",separators = "\t")

let ontology : Frame<string,string> = 
    let path = __SOURCE_DIRECTORY__ + @"..\..\data\Arabidopsis_Ontology.tab"
    Frame.ReadCsv(path,indexCol = "Key",separators = "\t")

let timeSeriesMatrix =
    timeSeriesData
    |> Frame.toArray2D
    |> Matrix.ofArray2D
    |> Matrix.transpose

let correlationMatrix = 
    Correlation.Matrix.columnWiseCorrelationMatrix Correlation.Seq.pearson timeSeriesMatrix
    |> Matrix.toArray2D

let (thr,stats) = FSharp.Stats.Testing.RMT.compute 0.9 0.001 0.01 correlationMatrix

let thresholdedMatrix = applyThreshold thr correlationMatrix

let thresholdedMatrixFrame =        
    thresholdedMatrix
    |> Frame.ofArray2D
    |> Frame.indexRowsWith timeSeriesData.RowKeys
    |> Frame.indexColsWith timeSeriesData.RowKeys

let finalNetwork = Frame.join JoinKind.Left thresholdedMatrixFrame ontology

let nodes : Series<string,string> =
    Frame.getCol "MapManNumber"  finalNetwork

let nodeConverter nodeLabel =
    [
        Grammar.Attribute.Label (nodeLabel);
    ]

Series.map (fun at mm -> Streamer.addNode nodeConverter at mm) nodes

let edgeConverter _ = 

    [
        Grammar.Attribute.EdgeType Grammar.EdgeDirection.Undirected
    ]

finalNetwork
|> Frame.map (fun prot1 prot2 corr -> 
    if prot1 <> prot2 && corr <> 0. then
        Streamer.addEdge edgeConverter (prot1+prot2) prot1 prot2 corr |> ignore
    )