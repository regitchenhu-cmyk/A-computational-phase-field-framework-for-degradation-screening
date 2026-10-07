$ErrorActionPreference = 'Stop'

$workspace = Split-Path -Parent $MyInvocation.MyCommand.Path
$matlab = (Get-Command matlab -ErrorAction Stop).Source
$logRoot = Join-Path $workspace 'Logs_Paper_ScientificallyCorrected_20260823\parallel_validation'
New-Item -ItemType Directory -Path $logRoot -Force | Out-Null

$cases = @('G50x70', 'dt20s', 'dt45s')
$jobs = foreach ($caseKey in $cases) {
    $outLog = Join-Path $logRoot ($caseKey + '.out.log')
    $errLog = Join-Path $logRoot ($caseKey + '.err.log')
    $batch = "cd('$($workspace.Replace("'", "''"))');run_one_scientifically_corrected_validation_case_20260823('$caseKey')"
    $suffix = switch ($caseKey) {
        'G50x70' { 'grid_50x70' }
        'dt20s' { 'dt_20s' }
        'dt45s' { 'dt_45s' }
        default { throw "Unknown validation case: $caseKey" }
    }
    $endFile = Join-Path $workspace ("Results_Paper_ScientificallyCorrected_20260823\Validation\cases\Standard_21MPa_80C_" + $suffix + "\end.mat")
    $process = Start-Process -FilePath $matlab -ArgumentList '-batch', $batch `
        -RedirectStandardOutput $outLog -RedirectStandardError $errLog `
        -WindowStyle Hidden -PassThru
    [pscustomobject]@{
        Case = $caseKey
        Process = $process
        Output = $outLog
        Error = $errLog
        EndFile = $endFile
    }
}

foreach ($job in $jobs) {
    $job.Process.WaitForExit()
    if ($job.Process.ExitCode -ne 0) {
        throw "Validation case $($job.Case) failed with exit code $($job.Process.ExitCode). See $($job.Error)"
    }
    if (-not (Test-Path -LiteralPath $job.EndFile -PathType Leaf)) {
        throw "Validation case $($job.Case) exited without the expected result: $($job.EndFile)"
    }
    Write-Output "DONE $($job.Case) exit=0"
}
