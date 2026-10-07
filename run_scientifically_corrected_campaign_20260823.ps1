$ErrorActionPreference = 'Stop'

$workspace = [System.IO.Path]::GetFullPath($PSScriptRoot)
$resultRoot = Join-Path $workspace 'Results_Parametric_ScientificallyCorrected_20260823'
$logRoot = Join-Path $workspace 'Logs_Parametric_ScientificallyCorrected_20260823'
$matlab = (Get-Command matlab -ErrorAction Stop).Source
$maxConcurrent = 4

$cases = @(
    'Grid_G18x26',
    'Grid_G28x40',
    'Grid_G40x56',
    'LengthScale_lh2.0',
    'LengthScale_lh3.0',
    'LengthScale_lh4.0',
    'LengthScale_lh5.0',
    'Orthogonal_O01_P21_T080_f020',
    'Orthogonal_O02_P21_T100_f050',
    'Orthogonal_O03_P21_T120_f100',
    'Orthogonal_O04_P28_T080_f050',
    'Orthogonal_O05_P28_T100_f100',
    'Orthogonal_O06_P28_T120_f020',
    'Orthogonal_O07_P35_T080_f100',
    'Orthogonal_O08_P35_T100_f020',
    'Orthogonal_O09_P35_T120_f050'
)

New-Item -ItemType Directory -Path $resultRoot -Force | Out-Null
New-Item -ItemType Directory -Path $logRoot -Force | Out-Null

$pending = [System.Collections.Generic.Queue[string]]::new()
foreach ($caseName in $cases) {
    $endFile = Join-Path (Join-Path $resultRoot $caseName) 'end.mat'
    if (Test-Path -LiteralPath $endFile) {
        Write-Output "SKIP completed $caseName"
    } else {
        $pending.Enqueue($caseName)
    }
}

$running = @()
$failures = @()
while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    while ($pending.Count -gt 0 -and $running.Count -lt $maxConcurrent) {
        $caseName = $pending.Dequeue()
        $stdout = Join-Path $logRoot ($caseName + '.out.log')
        $stderr = Join-Path $logRoot ($caseName + '.err.log')
        # Keep the -batch expression free of spaces because Start-Process
        # joins ArgumentList entries into a Windows command line.
        $expression = "cd('$($workspace.Replace("'", "''"))');run_scientifically_corrected_case_20260823('$caseName');"
        $process = Start-Process -FilePath $matlab `
            -ArgumentList @('-wait', '-batch', $expression) `
            -RedirectStandardOutput $stdout `
            -RedirectStandardError $stderr `
            -WindowStyle Hidden `
            -PassThru
        $running += [pscustomobject]@{
            Case = $caseName
            Process = $process
            Started = Get-Date
            Stdout = $stdout
            Stderr = $stderr
        }
        Write-Output "START $caseName PID=$($process.Id)"
    }

    Start-Sleep -Seconds 5
    $stillRunning = @()
    foreach ($job in $running) {
        if ($job.Process.HasExited) {
            $elapsed = (Get-Date) - $job.Started
            $endFile = Join-Path (Join-Path $resultRoot $job.Case) 'end.mat'
            if ($job.Process.ExitCode -eq 0 -and (Test-Path -LiteralPath $endFile)) {
                Write-Output ("DONE {0} exit=0 elapsed={1:hh\:mm\:ss}" -f $job.Case, $elapsed)
            } elseif ($job.Process.ExitCode -eq 0) {
                Write-Output ("FAIL {0} exit=0 but end.mat is missing elapsed={1:hh\:mm\:ss}" -f $job.Case, $elapsed)
                $failures += $job.Case
            } else {
                Write-Output ("FAIL {0} exit={1} elapsed={2:hh\:mm\:ss}" -f $job.Case, $job.Process.ExitCode, $elapsed)
                $failures += $job.Case
            }
        } else {
            $stillRunning += $job
        }
    }
    $running = $stillRunning
}

if ($failures.Count -gt 0) {
    throw ('Failed cases: ' + ($failures -join ', '))
}

Write-Output 'ALL_CASES_COMPLETE'
