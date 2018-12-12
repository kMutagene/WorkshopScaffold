// Include CsbScaffold
#load "../../.env/CsbScaffold.fsx"
#nowarn "10001"

open Deedle
open FSharp.Stats
open FSharpGephiStreamer
#time

///Filter out values of a matrix which are lower than the given threshold
let applyThreshold thr m = Array2D.map (fun x -> if abs(x) > thr then x else 0.) m

///Aggregated and filtered ratios of the proteins at the given timepoint
let timeSeriesData : Frame<string,string> = 
    let path = __SOURCE_DIRECTORY__ + @"..\..\data\AggregatedProteinQuantTable.tab"
    Frame.ReadCsv(path,indexCol = "Key",separators = "\t")

///MapMan Ontology Information on the proteins
let ontology : Frame<string,string> = 
    let path = __SOURCE_DIRECTORY__ + @"..\..\data\Arabidopsis_Ontology.tab"
    Frame.ReadCsv(path,indexCol = "Key",separators = "\t")

///Time Series Data in form of a matrix
let timeSeriesMatrix =
    timeSeriesData
    |> Frame.toArray2D
    |> Matrix.ofArray2D
    |> Matrix.transpose

///Pearson correlation values for all time series with all time series (Adjacency Matrix)
let correlationMatrix = 
    Correlation.Matrix.columnWiseCorrelationMatrix Correlation.Seq.pearson timeSeriesMatrix
    |> Matrix.toArray2D

///Critical Threshold computed by usage of RMT
let (thr,stats) = FSharp.Stats.Testing.RMT.compute 0.9 0.001 0.01 correlationMatrix

///Filtered Adjacency matrix 
let thresholdedMatrix = applyThreshold thr correlationMatrix

///Filtered Adjacency matrix in deedle data frame
let thresholdedMatrixFrame =        
    thresholdedMatrix
    |> Frame.ofArray2D
    |> Frame.indexRowsWith timeSeriesData.RowKeys
    |> Frame.indexColsWith timeSeriesData.RowKeys

///Final Result: Filtered Adjacency matrix in deedle data frame with appended ontology
let finalNetwork = Frame.join JoinKind.Left thresholdedMatrixFrame ontology

///Protein names with MapMan Number
let nodes : Series<string,string> =
    Frame.getCol "MapManNumber"  finalNetwork

///Gephi Converter for nodeLabel
let nodeConverter nodeLabel =
    [
        Grammar.Attribute.Label (nodeLabel);
    ]

//In this step we feed the nodes into gephi. Turn on streaming in gephi!
Series.map (fun at mm -> Streamer.addNode nodeConverter at mm) nodes

///Gephi converter for labeling edges as undirected
let edgeConverter _ = 

    [
        Grammar.Attribute.EdgeType Grammar.EdgeDirection.Undirected
    ]

//In this step we feed the edges into gephi.
finalNetwork
|> Frame.map (fun prot1 prot2 corr -> 
    if prot1 <> prot2 && corr <> 0. then
        Streamer.addEdge edgeConverter (prot1+prot2) prot1 prot2 corr |> ignore
    )