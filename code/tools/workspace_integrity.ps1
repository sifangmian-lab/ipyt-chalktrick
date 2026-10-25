# Workspace integrity-label repair for DSH Windows sandbox (with rollback).
#
# Symptom: confined (workspace-write) pwsh runs at Low integrity, while
# D:\vibelab\iypt\ipyt-chalktrick and its files were created at Medium integrity,
# so Windows no-write-up blocks every write there. The workspace root is already
# Low-labelled by DSH, but pre-existing subdirectories were not relabelled.
#
# This script only changes the mandatory integrity label (SACL). It does not
# change file contents, owners or DACLs. Modes:
#   Status    read-only report: token integrity, sample labels, object count
#   Pilot     temp-dir rehearsal: Low -> Medium -> delete (validates rollback)
#   Apply     Pilot, backup, recursive relabel to Low, verify, auto-rollback on failure
#   Rollback  recursive relabel back to Medium, verify
#
# Must run unconfined (full access); refuses a Low-integrity token.
# Administrator rights are NOT required; if Pilot's Low->Medium direction fails,
# do not run Apply and use an elevated PowerShell for rollback instead.

param(
    [ValidateSet('Status', 'Pilot', 'Apply', 'Rollback')]
    [string]$Mode = 'Status',
    [string]$Scope = 'D:\vibelab\iypt\ipyt-chalktrick',
    [string]$RecoveryDir = 'D:\vibelab\iypt\.integrity-recovery'
)

$ErrorActionPreference = 'Stop'

function Quote-Arg {
    param([string]$Value)
    if ($Value -match '\s') { return '"' + $Value + '"' }
    return $Value
}

function Invoke-Icacls {
    param([string[]]$Arguments)
    $id = [guid]::NewGuid().ToString('N')
    $outFile = Join-Path $env:TEMP "icacls-$id.out"
    $errFile = Join-Path $env:TEMP "icacls-$id.err"
    try {
        $proc = Start-Process -FilePath 'icacls.exe' `
            -ArgumentList (($Arguments | ForEach-Object { Quote-Arg $_ }) -join ' ') `
            -NoNewWindow -Wait -PassThru `
            -RedirectStandardOutput $outFile -RedirectStandardError $errFile
        $out = if (Test-Path -LiteralPath $outFile) { Get-Content -LiteralPath $outFile -Raw } else { '' }
        $err = if (Test-Path -LiteralPath $errFile) { Get-Content -LiteralPath $errFile -Raw } else { '' }
        return [pscustomobject]@{ ExitCode = $proc.ExitCode; Output = "$out`n$err" }
    }
    finally {
        Remove-Item -LiteralPath $outFile, $errFile -Force -ErrorAction SilentlyContinue
    }
}

function Get-IntegrityLabel {
    param([string]$Path)
    $res = Invoke-Icacls -Arguments @($Path)
    $m = [regex]::Match($res.Output, 'Mandatory Label\\(Low|Medium|High|System) Mandatory Level')
    if ($m.Success) { return $m.Groups[1].Value }
    return 'DEFAULT(=Medium)'
}

function Get-TokenIntegrity {
    $g = (whoami /groups 2>$null | Out-String)
    foreach ($lv in @('Low', 'Medium', 'High', 'System')) {
        if ($g -match "Mandatory Label\\$lv Mandatory Level") { return $lv }
    }
    return 'Unknown'
}

function Invoke-Pilot {
    $base = Split-Path -Parent $Scope
    if (-not (Test-Path -LiteralPath $base)) { throw "Pilot base directory not found: $base" }
    $dir = Join-Path $base ('.integrity-pilot-' + [guid]::NewGuid().ToString('N'))
    $sub = Join-Path $dir 'sub'
    New-Item -ItemType Directory -Path $sub -Force | Out-Null
    $file = Join-Path $sub 'pilot.txt'
    Set-Content -LiteralPath $file -Value 'integrity pilot' -Encoding UTF8
    try {
        $before = Get-IntegrityLabel $dir
        $r1 = Invoke-Icacls -Arguments @($dir, '/setintegritylevel', '(OI)(CI)L', '/T', '/C', '/Q')
        $dirLow = Get-IntegrityLabel $dir
        $fileLow = Get-IntegrityLabel $file
        $r2 = Invoke-Icacls -Arguments @($dir, '/setintegritylevel', '(OI)(CI)M', '/T', '/C', '/Q')
        $dirMed = Get-IntegrityLabel $dir
        $fileMed = Get-IntegrityLabel $file
        Write-Output "PILOT dir: $before -> L:$dirLow (exit=$($r1.ExitCode)) -> M:$dirMed (exit=$($r2.ExitCode))"
        Write-Output "PILOT file: L:$fileLow M:$fileMed"
        if ($r1.ExitCode -ne 0 -or $r2.ExitCode -ne 0) {
            throw "Pilot: icacls returned non-zero ($($r1.ExitCode)/$($r2.ExitCode))"
        }
        if ($dirLow -notmatch 'Low') { throw "Pilot: directory did not become Low" }
        if ($dirMed -notmatch 'Medium') {
            throw "Pilot: rollback direction (Low -> Medium) failed; this account cannot raise labels. Do NOT run Apply; use an elevated PowerShell for rollback."
        }
        if ($fileLow -notmatch 'Low' -or $fileMed -notmatch 'Medium') {
            throw "Pilot: file recursion/rollback failed (L:$fileLow M:$fileMed)"
        }
        Write-Output 'PILOT=OK'
    }
    finally {
        Remove-Item -LiteralPath $dir -Recurse -Force -ErrorAction SilentlyContinue
    }
}

function Get-SamplePaths {
    $candidates = @(
        $Scope,
        (Join-Path $Scope 'code\v2_model.py'),
        (Join-Path $Scope '.git\HEAD'),
        (Join-Path $Scope 'README.md')
    )
    return @($candidates | Where-Object { Test-Path -LiteralPath $_ })
}

function Invoke-Status {
    Write-Output "TOKEN_INTEGRITY=$(Get-TokenIntegrity)"
    Write-Output "SCOPE=$Scope EXISTS=$(Test-Path -LiteralPath $Scope)"
    if (Test-Path -LiteralPath $Scope) {
        $count = (Get-ChildItem -LiteralPath $Scope -Recurse -Force -ErrorAction SilentlyContinue |
                  Measure-Object).Count
        Write-Output "OBJECTS=$count"
        foreach ($s in Get-SamplePaths) { Write-Output "LABEL $s => $(Get-IntegrityLabel $s)" }
    }
    Write-Output "RECOVERYDIR=$RecoveryDir EXISTS=$(Test-Path -LiteralPath $RecoveryDir)"
    if (Test-Path -LiteralPath $RecoveryDir) {
        Get-ChildItem -LiteralPath $RecoveryDir -ErrorAction SilentlyContinue |
            Select-Object Name, Length, LastWriteTime |
            Format-Table -AutoSize | Out-String | Write-Output
    }
}

function Invoke-Apply {
    if (-not (Test-Path -LiteralPath $Scope)) { throw "Scope not found: $Scope" }
    $tok = Get-TokenIntegrity
    Write-Output "APPLY token=$tok scope=$Scope"
    if ($tok -eq 'Low') { throw "Current process is Low integrity; run unconfined (full access)." }
    Invoke-Pilot

    New-Item -ItemType Directory -Force -Path $RecoveryDir | Out-Null
    $stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
    $log = Join-Path $RecoveryDir "integrity-apply-$stamp.log"
    $rollback = Join-Path $RecoveryDir "ROLLBACK-$stamp.ps1"

    $pre = Invoke-Icacls -Arguments @($Scope)
    @(
        "APPLY START $(Get-Date -Format o)",
        "Scope=$Scope",
        "TokenIntegrity=$tok",
        "RecoveryDir=$RecoveryDir",
        "---- pre icacls ----",
        $pre.Output
    ) | Set-Content -LiteralPath $log -Encoding UTF8

    $rb = @"
# Auto-generated rollback: restore Medium integrity label on the subtree.
# Usage: powershell -NoProfile -ExecutionPolicy Bypass -File '$rollback'
`$ErrorActionPreference = 'Continue'
& icacls.exe '$Scope' '/setintegritylevel' '(OI)(CI)M' '/T' '/C' '/Q'
Write-Output "ROLLBACK_DONE exit=`$LASTEXITCODE"
"@
    Set-Content -LiteralPath $rollback -Value $rb -Encoding UTF8
    Write-Output "ROLLBACK_SCRIPT=$rollback"

    $apply = Invoke-Icacls -Arguments @($Scope, '/setintegritylevel', '(OI)(CI)L', '/T', '/C', '/Q')
    Add-Content -LiteralPath $log -Value @("---- apply icacls ----", $apply.Output)
    $failed = 0
    if ($apply.Output -match 'Failed processing (\d+) files') { $failed = [int]$Matches[1] }
    Write-Output "APPLY exit=$($apply.ExitCode) failed=$failed"

    $bad = @()
    foreach ($s in Get-SamplePaths) {
        $lab = Get-IntegrityLabel $s
        Write-Output "CHECK $s => $lab"
        Add-Content -LiteralPath $log -Value "CHECK $s => $lab"
        if ($lab -notmatch 'Low') { $bad += "$s => $lab" }
    }

    if ($apply.ExitCode -ne 0 -or $failed -gt 0 -or $bad.Count -gt 0) {
        Write-Output 'APPLY_FAILED -> AUTO ROLLBACK'
        $rbRes = Invoke-Icacls -Arguments @($Scope, '/setintegritylevel', '(OI)(CI)M', '/T', '/C', '/Q')
        Write-Output "AUTO_ROLLBACK exit=$($rbRes.ExitCode)"
        Add-Content -LiteralPath $log -Value @("---- auto rollback ----", $rbRes.Output)
        throw "Apply failed (exit=$($apply.ExitCode), failed=$failed, bad=$($bad -join '; ')); auto-rollback attempted. Manual rollback: powershell -NoProfile -File '$rollback'"
    }

    Write-Output 'APPLY=OK'
    Write-Output "LOG=$log"
    Write-Output "ROLLBACK_COMMAND=powershell -NoProfile -ExecutionPolicy Bypass -File '$rollback'"
}

function Invoke-Rollback {
    if (-not (Test-Path -LiteralPath $Scope)) { throw "Scope not found: $Scope" }
    $res = Invoke-Icacls -Arguments @($Scope, '/setintegritylevel', '(OI)(CI)M', '/T', '/C', '/Q')
    $failed = 0
    if ($res.Output -match 'Failed processing (\d+) files') { $failed = [int]$Matches[1] }
    Write-Output "ROLLBACK exit=$($res.ExitCode) failed=$failed"
    $bad = @()
    foreach ($s in Get-SamplePaths) {
        $lab = Get-IntegrityLabel $s
        Write-Output "CHECK $s => $lab"
        if ($lab -match 'Low') { $bad += "$s => $lab" }
    }
    if ($res.ExitCode -ne 0 -or $failed -gt 0 -or $bad.Count -gt 0) {
        throw "Rollback incomplete (exit=$($res.ExitCode), failed=$failed, stillLow=$($bad -join '; ')); try an elevated PowerShell with the same command"
    }
    Write-Output 'ROLLBACK=OK'
}

try {
    switch ($Mode) {
        'Status'   { Invoke-Status }
        'Pilot'    { Invoke-Pilot }
        'Apply'    { Invoke-Apply }
        'Rollback' { Invoke-Rollback }
    }
}
catch {
    Write-Output "ERROR: $($_.Exception.Message)"
    exit 1
}
exit 0
