using MultiSpinCoding
using Documenter

DocMeta.setdocmeta!(MultiSpinCoding, :DocTestSetup, :(using MultiSpinCoding); recursive=true)

makedocs(;
    modules=[MultiSpinCoding],
    authors="ArrogantGao <xz.gao@connect.ust.hk> and contributors",
    sitename="MultiSpinCoding.jl",
    format=Documenter.HTML(;
        canonical="https://ArrogantGao.github.io/MultiSpinCoding.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ArrogantGao/MultiSpinCoding.jl",
    devbranch="main",
)
