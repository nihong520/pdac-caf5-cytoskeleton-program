param(
    [string]$ProjectRoot = (Split-Path -Parent $PSScriptRoot),
    [switch]$Force
)

$ErrorActionPreference = "Stop"
$outDir = Join-Path $ProjectRoot "data\raw\pdacR"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$files = @(
    @{
        Name = "Moffitt_GEO_array.rds"
        Url = "https://raw.githubusercontent.com/rmoffitt/pdacR/master/data/Moffitt_GEO_array.rds"
        ExpectedMinimum = 17000000
    },
    @{
        Name = "TCGA_PAAD.rds"
        Url = "https://raw.githubusercontent.com/rmoffitt/pdacR/master/data/TCGA_PAAD.rds"
        ExpectedMinimum = 12000000
    },
    @{
        Name = "Puleo_array.rds"
        Url = "https://raw.githubusercontent.com/rmoffitt/pdacR/master/data/Puleo_array.rds"
        ExpectedMinimum = 41000000
    }
)

$records = @("file`turl`tsize_bytes`tsha256")
foreach ($item in $files) {
    $dest = Join-Path $outDir $item.Name
    $needsDownload = $Force -or !(Test-Path $dest)
    if (!$needsDownload -and (Get-Item $dest).Length -lt $item.ExpectedMinimum) {
        $needsDownload = $true
    }

    if ($needsDownload) {
        Write-Host "Downloading $($item.Name)..."
        & curl.exe --location --fail --retry 5 --retry-all-errors `
            --continue-at - --output $dest $item.Url
        if ($LASTEXITCODE -ne 0) {
            throw "Failed to download $($item.Name)."
        }
    } else {
        Write-Host "Using existing $($item.Name)."
    }

    $f = Get-Item $dest
    if ($f.Length -lt $item.ExpectedMinimum) {
        throw "$($item.Name) is unexpectedly small: $($f.Length) bytes."
    }
    try {
        $sha = (Get-FileHash -LiteralPath $dest -Algorithm SHA256 -ErrorAction Stop).Hash
    } catch {
        # Antivirus/indexing can transiently lock a newly completed large file.
        # Size and RDS integrity are checked separately by the analysis script.
        Write-Warning "Could not hash $($item.Name) because it is temporarily locked."
        $sha = "TEMPORARILY_LOCKED"
    }
    $records += "$($item.Name)`t$($item.Url)`t$($f.Length)`t$sha"
}

$records | Set-Content -LiteralPath (Join-Path $outDir "download_manifest.tsv") -Encoding UTF8
Write-Host "pdacR cohorts ready: $outDir"
