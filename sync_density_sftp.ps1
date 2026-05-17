param(
    [string]$RemoteRoot = "HPC_USERNAME@HPC_HOST:/home/w481/SCRATCH/miguelcores/",
    [string]$LocalRoot = "$PSScriptRoot"
)

# Parse user@host and remote path from user@host:/path
$parts = $RemoteRoot -split ":", 2
if ($parts.Length -ne 2) {
    Write-Error "RemoteRoot must be in the form user@host:/absolute/path"
    exit 1
}

$remoteHost = $parts[0]
$remotePath = $parts[1]

# List of concentration folders to sync
$concentrations = @(2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100)

# Build SFTP batch commands
$batchCommands = @()
$batchCommands += "lcd `"$LocalRoot`""
$batchCommands += "cd `"$remotePath`""

foreach ($conc in $concentrations) {
    $folderName = "${conc}_mmc"
    $localDir = Join-Path $LocalRoot "concentration_study\$folderName"
    
    # Ensure local directory exists
    if (!(Test-Path $localDir)) {
        New-Item -ItemType Directory -Path $localDir -Force | Out-Null
    }
    
    # Add command to pull density.xvg into the correct folder
    $batchCommands += "lcd `"$localDir`""
    $batchCommands += "cd `"$remotePath/concentration_study/$folderName`""
    $batchCommands += "get density.xvg"
}

$batchCommands += "bye"
$batch = $batchCommands -join "`n"

$tempBatch = New-TemporaryFile
Set-Content -Path $tempBatch -Value $batch -Encoding ASCII

Write-Host "Running sftp from $remoteHost ..."
Write-Host "Syncing density.xvg files for all concentration folders..."
& sftp -o StrictHostKeyChecking=accept-new -o PreferredAuthentications=password,keyboard-interactive -o PubkeyAuthentication=no -b $tempBatch $remoteHost

Remove-Item $tempBatch -Force

# Count how many files were pulled
$pulledCount = 0
$missingCount = 0
foreach ($conc in $concentrations) {
    $densityFile = Join-Path $LocalRoot "concentration_study\${conc}_mmc\density.xvg"
    if (Test-Path $densityFile) {
        $pulledCount++
    } else {
        $missingCount++
    }
}

Write-Host "`nSync complete: $pulledCount density.xvg files found, $missingCount missing"
