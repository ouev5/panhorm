var sankeyGoResults = [];
var sankeyKeggResults = [];
var currentSankeySource = 'kegg';

function runSankeyEnrich(){
    if(!currentElements || currentElements.length === 0){
        setStatus('Please select genes first','error');
        return;
    }
    setStatus('Running GO and KEGG enrichment...','success');
    document.getElementById('sankeyCard').style.display = 'block';
    document.getElementById('sankeyCard').scrollIntoView({behavior:'smooth'});
    
    // 同时执行 GO 和 KEGG 富集
    var goPromise = fetch('/api/enrichment', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({genes: currentElements, type: 'go'})
    }).then(r => r.json());
    
    var keggPromise = fetch('/api/enrichment', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({genes: currentElements, type: 'kegg'})
    }).then(r => r.json());
    
    Promise.all([goPromise, keggPromise])
    .then(function(data){
        var goData = data[0];
        var keggData = data[1];
        
        sankeyGoResults = goData.results || [];
        sankeyKeggResults = keggData.results || [];
        
        // 默认显示 KEGG
        currentSankeySource = 'kegg';
        
        if(sankeyKeggResults.length > 0 || sankeyGoResults.length > 0){
            renderSankeyDot(sankeyKeggResults, 'kegg');
            setStatus('GO: '+sankeyGoResults.length+' pathways, KEGG: '+sankeyKeggResults.length+' pathways','success');
        } else {
            setStatus('No pathways found','error');
        }
    })
    .catch(function(e){
        setStatus('Error: '+e.message,'error');
    });
}

function showContinueModal(){
    var results = currentSankeySource === 'go' ? sankeyGoResults : sankeyKeggResults;
    var allGenes = [];
    if(results && results.length > 0){
        results.forEach(function(r){
            if(r.genes){
                var genes = Array.isArray(r.genes) ? r.genes : r.genes.split(';');
                genes.forEach(function(g){
                    if(g && !allGenes.includes(g)) allGenes.push(g);
                });
            }
        });
    }
    if(allGenes.length === 0){setStatus('No genes found','error');return;}
    
    var geneCheckboxes = allGenes.slice(0, 50).map(function(g, i){
        return '<label style="display:inline-flex;align-items:center;gap:5px;padding:4px 8px;margin:2px;background:#f1f5f9;border-radius:4px;font-size:0.8rem;cursor:pointer;">'+
            '<input type="checkbox" checked onchange="toggleGeneSel(this,\''+g+'\')"> '+g+'</label>';
    }).join('');
    
    var html = '<div id="cm" style="position:fixed;z-index:9999;left:0;top:0;width:100%;height:100%;background:rgba(0,0,0,0.5);display:flex;align-items:center;justify-content:center;">'+
        '<div style="background:white;padding:20px;border-radius:12px;max-width:600px;max-height:80vh;overflow-y:auto;">'+
        '<h3 style="margin:0 0 10px;"><i class="fa fa-check-square"></i> Select Genes for Next Step</h3>'+
        '<p style="color:#64748b;font-size:0.85rem;margin-bottom:10px;">Found '+allGenes.length+' genes from enrichment. Select genes to continue:</p>'+
        '<div style="margin-bottom:10px;">'+
        '<button class="btn btn-secondary" onclick="selectAllGenes(true)" style="font-size:0.7rem;padding:2px 8px;margin-right:5px;">Select All</button>'+
        '<button class="btn btn-secondary" onclick="selectAllGenes(false)" style="font-size:0.7rem;padding:2px 8px;">Deselect All</button>'+
        '</div>'+
        '<div id="geneSelBox" style="margin-bottom:15px;padding:10px;background:#f8fafc;border-radius:8px;max-height:250px;overflow-y:auto;">'+geneCheckboxes+'</div>'+
        '<div style="display:flex;gap:8px;">'+
        '<button class="btn btn-accent" onclick="goToStringSel()" style="flex:1;"><i class="fa fa-share-alt"></i> STRING Database</button>'+
        '<button class="btn btn-primary" onclick="showSCModalSel()" style="flex:1;"><i class="fa fa-cells"></i> Single-Cell Map</button>'+
        '</div>'+
        '<button class="btn btn-secondary" onclick="document.getElementById(\'cm\').remove()" style="width:100%;margin-top:10px;">Cancel</button>'+
        '</div></div>';
    
    selectedSankeyGenes = allGenes.slice(0, 50);
    document.body.insertAdjacentHTML('beforeend', html);
}

var selectedSankeyGenes = [];

function toggleGeneSel(cb, gene){
    if(cb.checked){if(!selectedSankeyGenes.includes(gene))selectedSankeyGenes.push(gene);}
    else{var i=selectedSankeyGenes.indexOf(gene);if(i>-1)selectedSankeyGenes.splice(i,1);}
}

function selectAllGenes(selAll){
    var cbs = document.getElementById('geneSelBox').querySelectorAll('input[type="checkbox"]');
    cbs.forEach(function(cb){
        cb.checked = selAll;
        var gene = cb.parentElement.textContent.trim();
        if(selAll && !selectedSankeyGenes.includes(gene)) selectedSankeyGenes.push(gene);
        else if(!selAll) selectedSankeyGenes = [];
    });
}

function goToStringSel(){
    if(selectedSankeyGenes.length===0){alert('Please select at least one gene');return;}
    fetch('/api/string_url', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({genes: selectedSankeyGenes})
    })
    .then(function(r){return r.json()})
    .then(function(d){
        if(d.url){
            document.getElementById('cm').remove();
            window.open(d.url, '_blank');
        } else {
            alert('Failed to generate STRING URL');
        }
    })
    .catch(function(e){
        alert('Error: ' + e.message);
    });
}

function showSCModalSel(){
    if(selectedSankeyGenes.length===0){alert('Please select at least one gene');return;}
    document.getElementById('cm').remove();
    var html='<div id="scm" style="position:fixed;z-index:9999;left:0;top:0;width:100%;height:100%;background:rgba(0,0,0,0.5);display:flex;align-items:center;justify-content:center;">'+
        '<div style="background:white;padding:20px;border-radius:12px;max-width:450px;">'+
        '<h3 style="margin:0 0 15px;"><i class="fa fa-cells"></i> Single-Cell Map Options</h3>'+
        '<p style="color:#64748b;font-size:0.85rem;margin-bottom:15px;">'+selectedSankeyGenes.length+' genes selected. Choose data source:</p>'+
        '<button class="btn btn-primary" onclick="uploadOwnMatrix()" style="width:100%;margin-bottom:8px;text-align:left;padding:12px;"><i class="fa fa-upload"></i> <b>Upload Your Matrix</b><br><span style="font-size:0.75rem;color:#94a3b8;">Upload your own single-cell expression matrix</span></button>'+
        '<button class="btn btn-secondary" onclick="openPublicDB()" style="width:100%;margin-bottom:8px;text-align:left;padding:12px;"><i class="fa fa-database"></i> <b>Public Databases</b><br><span style="font-size:0.75rem;color:#94a3b8;">GEO, EMBL-EBI, Human Cell Atlas</span></button>'+
        '<button class="btn btn-success" onclick="useLocalData()" style="width:100%;margin-bottom:8px;text-align:left;padding:12px;"><i class="fa fa-check-circle"></i> <b>Use Local Data</b> <span style="background:#10b981;color:white;padding:2px 6px;border-radius:4px;font-size:0.65rem;">Tabula Sapiens</span><br><span style="font-size:0.75rem;color:#94a3b8;">94,836 cells, 24 tissues (pre-loaded)</span></button>'+
        '<button class="btn btn-secondary" onclick="document.getElementById(\'scm\').remove()" style="width:100%;margin-top:10px;">Cancel</button></div></div>';
    document.body.insertAdjacentHTML('beforeend', html);
}

function uploadOwnMatrix(){
    alert('Matrix upload feature coming soon!\\nSupported formats: .h5ad, .csv, .mtx');
}

function openPublicDB(){
    var html='<div id="pdb" style="position:fixed;z-index:9999;left:0;top:0;width:100%;height:100%;background:rgba(0,0,0,0.5);display:flex;align-items:center;justify-content:center;">'+
        '<div style="background:white;padding:20px;border-radius:12px;max-width:500px;">'+
        '<h3 style="margin:0 0 15px;"><i class="fa fa-external-link"></i> Public Single-Cell Databases</h3>'+
        '<a href="https://www.ncbi.nlm.nih.gov/geo/" target="_blank" style="display:block;padding:10px;border:1px solid #e2e8f0;border-radius:8px;margin-bottom:8px;text-decoration:none;color:inherit;"><b>GEO (Gene Expression Omnibus)</b><br><span style="font-size:0.75rem;color:#64748b;">NCBI public repository</span></a>'+
        '<a href="https://www.ebi.ac.uk/arrayexpress/" target="_blank" style="display:block;padding:10px;border:1px solid #e2e8f0;border-radius:8px;margin-bottom:8px;text-decoration:none;color:inherit;"><b>EMBL-EBI ArrayExpress</b><br><span style="font-size:0.75rem;color:#64748b;">European bioinformatics institute</span></a>'+
        '<a href="https://www.humancellatlas.org/" target="_blank" style="display:block;padding:10px;border:1px solid #e2e8f0;border-radius:8px;margin-bottom:8px;text-decoration:none;color:inherit;"><b>Human Cell Atlas</b><br><span style="font-size:0.75rem;color:#64748b;">Comprehensive cell reference map</span></a>'+
        '<a href="https://tabula-sapiens.sf.czbiohub.org/" target="_blank" style="display:block;padding:10px;border:1px solid #e2e8f0;border-radius:8px;margin-bottom:8px;text-decoration:none;color:inherit;"><b>Tabula Sapiens</b><br><span style="font-size:0.75rem;color:#64748b;">Human cell atlas from CZ Biohub</span></a>'+
        '<button class="btn btn-secondary" onclick="document.getElementById(\'pdb\').remove()" style="width:100%;margin-top:10px;">Close</button></div></div>';
    document.body.insertAdjacentHTML('beforeend', html);
}

function useLocalData(){
    document.getElementById('scm').remove();
    currentElements = selectedSankeyGenes;
    showSingleCellCluster();
}
