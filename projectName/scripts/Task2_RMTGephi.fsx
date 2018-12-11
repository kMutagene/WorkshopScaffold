// Include CsbScaffold
#load "../../.env/CsbScaffold.fsx"
#load "Task1_Deedle.fsx"
// If you want to use the wrappers for unmanaged LAPACK functions from of FSharp.Stats 
// include the path to the .lib folder manually to your PATH environment variable and make sure you set FSI to 64 bit

// use the following lines of code to ensure that LAPACK functionalities are enabled if you want to use them
// fails with "MKL service either not available, or not started" if lib folder is not included in PATH.
//open FSharp.Stats
//FSharp.Stats.Algebra.LinearAlgebra.Service()
#nowarn 
open System
open FSharpAux
open Task1_Deedle
open Deedle
open FSharp.Stats
open FSharpGephiStreamer
#time
let y : Frame<string,string> = 
    Frame.ReadCsv(@"E:\Users\Lukas\Source\Repos\WorkshopScaffold\projectName\data\AggregatedProteinQuantTable.tab",indexCol = "Key",separators = "\t")

let ontology = 
    let path = @"E:\Users\Lukas\OneDrive - tukl\CSB\CsbScaffold\Berlin Workshop\data\proteinQuants_Filtered_Normalized_SetNaN.tab"
    use sr = System.IO.StreamReader(path)
    sr.ReadToEnd()
    |> String.split '\n'
    |> Array.map (fun x -> 
        let cancer = String.split '\t' x
        let mapmanNumber = 
            if cancer.[1] = "" then 
                ""
            else 
                cancer.[1] |> String.skip 4 |> String.split '.' |> fun x -> 
                    if x.Length < 3 then "" 
                    else
                        x.[0] + "." + x.[1]  +  x.[2]
        cancer.[0],(mapmanNumber,cancer.[2])
        )
    |> Map.ofArray


let sameLabelCount = 

    ontology
    |> Map.toArray
    |> Array.countBy (fun (id,(md,_)) -> md)
    |> Array.filter (fun (x,y) -> x = "" |> not)
    |> Array.sortByDescending snd
    |> Array.take 15
    |> Array.map fst

let onlyBiggestGroups = 
    y
    |> Frame.filterRows (fun key _ -> 
        
        let b = Array.contains (fst ontology.[key]) sameLabelCount
        //printfn "ID: %s; mapman: %s; foundIT: %b" key (fst ontology.[key]) b
        b)

let matrix =
    onlyBiggestGroups
    |> Frame.toArray2D
    //|> Array2D.toJaggedArray
    //|> Array.take 300  
    //|> Array2D.ofJaggedArray
    |> Matrix.ofArray2D
    |> Matrix.transpose

meanAcrossBiologicalReplicates.ColumnKeys
|> Seq.toArray

let keys = onlyBiggestGroups.RowKeys |> Seq.toArray

let corr = 
    Correlation.Matrix.columnWiseCorrelationMatrix Correlation.Seq.pearson matrix
#time
//0.9873 

//let (thr,stats) = FSharp.Stats.Testing.RMT.compute 0.9 0.001 0.01 (Matrix.toArray2D corr)
let thr = 0.994140625

let thresholded = Matrix.map (fun x -> if x > thr then x else 0.) corr

for i = 0 to (matrix.Dimensions |> snd) - 1 do Streamer.addNodeBy string i |> ignore

let nodeConverter nodeId =
    [
        Grammar.Attribute.Label (fst ontology.[keys.[nodeId]]);
    ]

let edgeConverter edgeId = 
    [
        Grammar.Attribute.EdgeType Grammar.EdgeDirection.Undirected
    ]    

for i = 0 to (matrix.Dimensions |> snd) - 1 do
    Streamer.updateNode nodeConverter i i

let mutable edges = 0
for i = 0 to (thresholded.Dimensions |> fst) - 1 do
    for j = i + 1 to (thresholded.Dimensions |> fst) - 1 do
        let x = thresholded.[i,j]
        if x = 0. |> not then 
            Streamer.addEdge edgeConverter edges i j x |> ignore
            edges <- edges + 1    




matrix.Dimensions
corr.Dimensions