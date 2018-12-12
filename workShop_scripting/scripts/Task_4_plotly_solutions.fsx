// Include CsbScaffold
#nowarn "10001"
#load "../../.env/CsbScaffold.fsx"
#load "Task_3_2_Deedle.fsx"
// If you want to use the wrappers for unmanaged LAPACK functions from of FSharp.Stats 
// include the path to the .lib folder manually to your PATH environment variable and make sure you set FSI to 64 bit

// use the following lines of code to ensure that LAPACK functionalities are enabled if you want to use them
// fails with "MKL service either not available, or not started" if lib folder is not included in PATH.
//open FSharp.Stats
//FSharp.Stats.Algebra.LinearAlgebra.Service()

open BioFSharp
open Deedle
open FSharpAux
open FSharpAux.IO
open System.IO
open Task_3_2_Deedle
open FSharp.Stats
open FSharp.Plotly

//Task 1: For each timepoint, plot the distribution of the mean values for all proteins 
meanAcrossBiologicalReplicates
|> Frame.getNumericCols 
|> Series.map (fun colkey colSeries -> Chart.Histogram(colSeries |> Series.values, Name=colkey) )
|> Series.values
|> Chart.Stack (2,0.1)
|> Chart.withSize (1500.,1500.)
|> Chart.Show


//Task 2: For a proteins of your interest, create a range plot or line plot showing the time course of its abundance and its dispersion.
let meanPoi : float [] =
    Frame.getRow "AT1G02920.1" meanAcrossBiologicalReplicates
    |> Series.values
    |> Array.ofSeq

let stdevPoiPos : float [] =
    Frame.getRow "AT1G02920.1" stdevAcrossBiologicalReplicates
    |> Series.values
    |> Array.ofSeq
    |> Array.mapi (fun i x -> meanPoi.[i] + x)

let stdevPoiNeg : float [] =
    Frame.getRow "AT1G02920.1" stdevAcrossBiologicalReplicates
    |> Series.values
    |> Array.ofSeq
    |> Array.mapi (fun i x -> meanPoi.[i] - x)

Chart.Range([1;2;3;4;5],meanPoi,stdevPoiPos,stdevPoiNeg,Color="grey",RangeColor="lightblue",Name="AT1G02920.1")
|> Chart.withTitle("Time course of abundance and dispersion of AT1G02920.1")
|> Chart.withX_AxisStyle("TimePoint",Showline=true,Showgrid=false)
|> Chart.withY_AxisStyle("log2meanRatio(N14/N15)",Showline=true,Showgrid=false)
|> Chart.Show

