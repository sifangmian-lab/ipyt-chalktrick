[CmdletBinding()]
param(
    [string]$LocalPath = 'D:\vibelab\iypt',
    [string]$Branch = 'develop'
)

$ErrorActionPreference = 'Stop'
$Repository = 'https://github.com/sifangmian-lab/ipyt-chalktrick.git'

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw '未检测到 Git。请先安装 Git for Windows，再重新运行本脚本。'
}

$FullPath = [System.IO.Path]::GetFullPath($LocalPath)
$GitDir = Join-Path $FullPath '.git'

if (-not (Test-Path -LiteralPath $FullPath)) {
    $Parent = Split-Path -Parent $FullPath
    New-Item -ItemType Directory -Force -Path $Parent | Out-Null

    & git clone --branch $Branch --single-branch $Repository $FullPath
    if ($LASTEXITCODE -ne 0) {
        throw "首次克隆失败；目标路径：$FullPath"
    }
}
elseif (-not (Test-Path -LiteralPath $GitDir)) {
    $ExistingItems = Get-ChildItem -LiteralPath $FullPath -Force
    if ($ExistingItems) {
        throw "目标目录已存在且不是空的 Git 工作区。为避免覆盖现有文件，脚本已停止：$FullPath"
    }

    & git clone --branch $Branch --single-branch $Repository $FullPath
    if ($LASTEXITCODE -ne 0) {
        throw "首次克隆失败；目标路径：$FullPath"
    }
}
else {
    $RemoteUrl = (& git -C $FullPath remote get-url origin 2>$null)
    if ($LASTEXITCODE -ne 0) {
        throw '现有工作区缺少 origin 远程仓库，脚本未做修改。'
    }
    if ($RemoteUrl -notmatch 'github\.com[:/]sifangmian-lab/ipyt-chalktrick(?:\.git)?$') {
        throw "现有工作区的 origin 不是本项目，脚本已停止：$RemoteUrl"
    }

    & git -C $FullPath fetch origin $Branch
    if ($LASTEXITCODE -ne 0) {
        throw '获取远端更新失败。'
    }

    & git -C $FullPath checkout $Branch
    if ($LASTEXITCODE -ne 0) {
        throw "无法切换到 $Branch 分支；请先处理本地未提交修改。"
    }

    & git -C $FullPath pull --ff-only origin $Branch
    if ($LASTEXITCODE -ne 0) {
        throw '无法快进同步。可能存在本地提交或分支分歧；脚本未强制覆盖，请检查 git status。'
    }
}

Write-Host ''
Write-Host "工作区路径：$FullPath"
Write-Host "跟踪分支：$Branch"
& git -C $FullPath status --short --branch
Write-Host ''
Write-Host '同步完成。脚本不会自动提交或推送本地修改。'
