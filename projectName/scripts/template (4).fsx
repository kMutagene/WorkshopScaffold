// Include CsbScaffold
#load "../../.env/CsbScaffold.fsx"
// If you want to use the wrappers for unmanaged LAPACK functions from of FSharp.Stats 
// include the path to the .lib folder manually to your PATH environment variable and make sure you set FSI to 64 bit

// use the following lines of code to ensure that LAPACK functionalities are enabled if you want to use them
// fails with "MKL service either not available, or not started" if lib folder is not included in PATH.
//open FSharp.Stats
//FSharp.Stats.Algebra.LinearAlgebra.Service()
#nowarn 
open BioFSharp
open Deedle
open FSharpAux
open FSharpAux.IO
open System.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute

type DataSchema = {

    Protein     : string
    [<FieldAttribute(0)>]
    MapManNumer : string
    [<FieldAttribute(1)>]
    Description : string
    [<FieldAttribute(2)>]
    NA_1_1      : float
    [<FieldAttribute(3)>]
    NA_1_2      : float
    [<FieldAttribute(4)>]
    NA_1_3      : float
    [<FieldAttribute(5)>]
    NA_2_1      : float
    [<FieldAttribute(6)>]
    NA_2_2      : float
    [<FieldAttribute(7)>]
    NA_2_3      : float
    [<FieldAttribute(8)>]
    NA_3_1      : float
    [<FieldAttribute(9)>]
    NA_3_2      : float
    [<FieldAttribute(10)>]
    NA_3_3      : float
    [<FieldAttribute(11)>]
    CA_1_1      : float
    [<FieldAttribute(12)>]
    CA_1_2      : float
    [<FieldAttribute(13)>]
    CA_1_3      : float
    [<FieldAttribute(14)>]
    CA_2_1      : float
    [<FieldAttribute(15)>]
    CA_2_2      : float
    [<FieldAttribute(16)>]
    CA_2_3      : float
    [<FieldAttribute(17)>]
    CA_3_1      : float
    [<FieldAttribute(18)>]
    CA_3_2      : float
    [<FieldAttribute(19)>]
    CA_3_3      : float
    [<FieldAttribute(20)>]
    DA6_1_1     : float
    [<FieldAttribute(21)>]
    DA6_1_2     : float
    [<FieldAttribute(22)>]
    DA6_1_3     : float
    [<FieldAttribute(23)>]
    DA6_2_1     : float
    [<FieldAttribute(24)>]
    DA6_2_2     : float
    [<FieldAttribute(25)>]
    DA6_2_3     : float
    [<FieldAttribute(26)>]
    DA6_3_1     : float
    [<FieldAttribute(27)>]
    DA6_3_2     : float
    [<FieldAttribute(28)>]
    DA6_3_3     : float
    [<FieldAttribute(29)>]
    DA12_1_1    : float
    [<FieldAttribute(30)>]
    DA12_1_2    : float
    [<FieldAttribute(31)>]
    DA12_1_3    : float
    [<FieldAttribute(32)>]
    DA12_2_1    : float
    [<FieldAttribute(33)>]
    DA12_2_2    : float
    [<FieldAttribute(34)>]
    DA12_2_3    : float
    [<FieldAttribute(35)>]
    DA12_3_1    : float
    [<FieldAttribute(36)>]
    DA12_3_2    : float
    [<FieldAttribute(37)>]
    DA12_3_3    : float
    [<FieldAttribute(38)>]
    DA24_1_1    : float
    [<FieldAttribute(39)>]
    DA24_1_2    : float
    [<FieldAttribute(40)>]
    DA24_1_3    : float
    [<FieldAttribute(41)>]
    DA24_2_1    : float
    [<FieldAttribute(42)>]
    DA24_2_2    : float
    [<FieldAttribute(43)>]
    DA24_2_3    : float
    [<FieldAttribute(44)>]
    DA24_3_1    : float
    [<FieldAttribute(45)>]
    DA24_3_2    : float
    [<FieldAttribute(46)>]
    DA24_3_3    : float
}

let filePaths = 
    @"C:\Users\Kevin\source\repos\WorkshopScaffold\projectName\data\ReplicateCorrs"
    |> Directory.EnumerateFiles
    |> Seq.filter (fun x -> (Path.GetExtension x) = ".tab")
    |> Array.ofSeq

let reader = new Csv.CsvReader<DataSchema>(schemaMode=Csv.Fill)

let rawData =
    filePaths
    |> Array.map (fun name -> name,reader.ReadFile(name,'\t',false,1) |> Array.ofSeq)
