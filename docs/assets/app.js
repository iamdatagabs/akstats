const DATA_URL = "data/dashboard.json";
let dashboard = null, mafChart = null, deltaChart = null, patientMafChart = null, prevalenceChart = null, geoExpressionChart = null, geoQcChart = null, chartMode = "bar";
const $ = (id) => document.getElementById(id);
const stageNames = {Normal:"Pele normal",AK:"Queratose actínica",SCC:"Carcinoma espinocelular",cSCC:"Carcinoma espinocelular"};
const fmt = (n,digits=2) => Number(n).toLocaleString("pt-BR",{maximumFractionDigits:digits});
const pct = (n) => `${fmt(Number(n)*100,1)}%`;
const pFmt = (n) => n===null || n===undefined || n==="" || !Number.isFinite(Number(n)) ? "não estimável" : (Number(n)<.001 ? "< 0,001" : fmt(n,3));
const escapeHtml = (s) => String(s??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
function showStatus(message,error=false){const el=$("data-status");el.textContent=message;el.classList.toggle("error",error);}
function fillKpis(s){$("kpi-genes").textContent=fmt(s.panel_genes??0,0);$("kpi-raw").textContent=fmt(s.source_variants??0,0);$("kpi-axis").textContent=fmt(s.axis_variants??0,0);$("kpi-patients").textContent=fmt(s.patients??0,0);const notes=["Número de genes contemplados no painel de sequenciamento; não representa genes alterados.","Registros de variantes na tabela suplementar antes do filtro funcional.","Variantes funcionais retidas no eixo ARID2–MAPK após a curadoria e o filtro do painel.","Indivíduos representados na tabela de origem; a unidade para as novas análises pareadas."];document.querySelectorAll(".kpi small").forEach((el,i)=>{if(notes[i])el.textContent=notes[i];});}
function filteredRows(){
 if(!dashboard?.maf)return [];
 const stage=$("stage-filter").value,module=$("module-filter").value,query=$("gene-search").value.trim().toLocaleLowerCase("pt-BR");
 return dashboard.maf.filter(r=>(stage==="all"||r.Stage===stage)&&(module==="all"||r.Module===module)&&(!query||[r.Module,r.Biological_group,r.Stage].join(" ").toLocaleLowerCase("pt-BR").includes(query)));
}
function drawMaf(){
 const rows=filteredRows(),stages=[...new Set(rows.map(r=>r.Stage))],labels=stages.map(s=>stageNames[s]||s);
 const datasets=["ARID2","MAPK"].map((group,i)=>({label:group,data:stages.map(s=>{const item=rows.find(r=>r.Stage===s&&r.Biological_group===group);return item?Number(item.median)*100:null;}),backgroundColor:i?"#879675":"#a67c59",borderColor:i?"#526b4b":"#76543d",borderWidth:2,borderRadius:6,tension:.3,spanGaps:true}));
 if(mafChart)mafChart.destroy();
 mafChart=new Chart($("maf-chart"),{type:chartMode,data:{labels,datasets},options:{responsive:true,interaction:{mode:"index",intersect:false},plugins:{legend:{position:"bottom"},tooltip:{callbacks:{label:c:`${c.dataset.label}: ${fmt(c.raw)}%`}}},scales:{y:{beginAtZero:true,title:{display:true,text:"MAF mediana (%)"},grid:{color:"#eee8de"}},x:{grid:{display:false}}},onClick:(event,elements)=>{if(!elements.length)return;$("stage-filter").value=stages[elements[0].index];renderExplorer();}}});
 addChartNote("maf-chart","MAF é a fração de leituras com a variante. A mediana resume os registros de variantes por grupo e estágio; não equivale a tamanho clonal nem a desfecho clínico.");
 $("result-count").textContent=`${rows.length} registro(s) exibido(s) · ${new Set(rows.map(r=>r.Module)).size} módulo(s)`;
 renderTable(rows);
}
function addChartNote(canvasId,text){const canvas=$(canvasId);if(!canvas)return;const card=canvas.closest(".chart-card");if(!card||card.querySelector(".metric-explanation"))return;const note=document.createElement("p");note.className="metric-explanation";note.textContent=text;card.append(note);}
function drawDelta(){
 const rows=dashboard?.modules||[];if(deltaChart)deltaChart.destroy();
 deltaChart=new Chart($("delta-chart"),{type:"bar",data:{labels:rows.map(r=>String(r.Module||"").replace(/^\\d+\\. /,"")),datasets:[{label:"ΔMAF mediana (p.p.)",data:rows.map(r=>Number(r.median_delta_maf)*100),backgroundColor:"#879675",borderRadius:5}]},options:{indexAxis:"y",responsive:true,plugins:{legend:{display:false},tooltip:{callbacks:{label:c:`${fmt(c.raw)} pontos percentuais`}}},scales:{x:{title:{display:true,text:"Pontos percentuais"},grid:{color:"#eee8de"}},y:{grid:{display:false}}}}});
}
function drawCooc(){
 addChartNote("delta-chart","ΔMAF = MAF em cSCC menos MAF em AK. Valores positivos indicam fração alélica observada maior em cSCC; negativos indicam menor. Não demonstram expansão clonal.");
 const host=$("cooc-summary");host.replaceChildren();
 (dashboard?.fisher||[]).forEach(r=>{const row=document.createElement("div");row.className="cooc-row";row.innerHTML=`<span>${escapeHtml(stageNames[r.stage]||r.stage)}</span><b>${fmt(r.both_present,0)} / ${fmt((r.both_present||0)+(r.arid2_only||0)+(r.mapk_only||0)+(r.neither||0),0)} pacientes</b>`;host.append(row);
 const note=document.createElement("p");note.className="metric-explanation";note.textContent=`Coocorrência observada: ${fmt(r.both_present,0)} pacientes com alterações ARID2 e MAPK acima de 2% no mesmo estágio. A tabela 2×2 e o teste exato de Fisher descrevem associação, não interação causal.`;host.append(note);});
 if(!host.children.length)host.textContent="Sem registros de coocorrência disponíveis.";
}
function renderTable(rows){
 const body=$("maf-table");body.replaceChildren();
 if(!rows.length){const tr=document.createElement("tr");tr.innerHTML='<td colspan="4">Nenhum registro corresponde aos filtros atuais.</td>';body.append(tr);return;}
 rows.slice().sort((a,b)=>Number(b.median)-Number(a.median)).forEach(r=>{const tr=document.createElement("tr");tr.innerHTML=`<td>${escapeHtml(stageNames[r.Stage]||r.Stage)}</td><td><span class="group-pill ${r.Biological_group==="ARID2"?"arid":"mapk"}">${escapeHtml(r.Biological_group)}</span></td><td>${escapeHtml(String(r.Module||"—").replace(/^\\d+\\. /,""))}</td><td class="numeric">${pct(r.median)}</td>`;body.append(tr);});
}
function ensurePatientSection(){
 if($("patient-analysis"))return;
 const section=document.createElement("section");section.id="patient-analysis";section.className="patient-analysis";
 section.innerHTML=`<div class="patient-heading"><div><span class="eyebrow">ANÁLISE POR UNIDADE BIOLÓGICA</span><h3>Métricas no nível do <em>paciente.</em></h3><p>Agregação por paciente evita tratar variantes repetidas do mesmo indivíduo como observações independentes.</p></div><span class="live-tag">n = pacientes</span></div>
 <div class="patient-chart-grid"><article class="chart-card"><div class="chart-title"><div><span>MAF AGREGADA</span><h3>Mediana da MAF máxima por paciente</h3></div></div><canvas id="patient-maf-chart" height="170"></canvas><p class="metric-explanation">Para cada paciente, grupo e estágio, considera-se a maior MAF entre as variantes funcionais observadas. A mediana resume a distribuição entre pacientes; ausência de variante elegível é codificada como zero observado no conjunto filtrado.</p></article>
 <article class="chart-card"><div class="chart-title"><div><span>PREVALÊNCIA MOLECULAR</span><h3>Pacientes acima do limiar de 2%</h3></div></div><canvas id="prevalence-chart" height="170"></canvas><p class="metric-explanation">Percentual de pacientes com ao menos uma variante do grupo com MAF ≥ 2%. As barras de erro representam IC 95% de Wilson, apropriado para proporções em amostras pequenas; o limiar é operacional, não clínico.</p></article></div>
 <article class="chart-card patient-tests"><div class="chart-title"><div><span>TESTES PAREADOS</span><h3>Variação entre AK e cSCC</h3></div></div><div class="table-scroll"><table class="data-table"><thead><tr><th>Grupo</th><th>n</th><th>Friedman p</th><th>Wilcoxon p</th><th>Ganho / perda</th><th>McNemar exato p</th></tr></thead><tbody id="patient-tests-body"></tbody></table></div><p class="metric-explanation">Friedman compara os três estágios de forma pareada; Wilcoxon compara AK e cSCC; McNemar exato avalia mudanças na presença acima de 2%. São análises exploratórias, com baixa precisão esperada em uma coorte de 16 pacientes. p-valor não mede magnitude nem relevância biológica.</p></article></section>`;
 $("maf-table").closest(".data-table-card").before(section);
}
function drawPatientMetrics(){
 ensurePatientSection();
 const rows=dashboard?.patient_metrics||[],groups=["ARID2","MAPK"],stages=["Normal","AK","cSCC"];
 const selected=(g,s)=>rows.find(r=>r.Biological_group===g&&r.Stage===s);
 if(patientMafChart)patientMafChart.destroy();
 patientMafChart=new Chart($("patient-maf-chart"),{type:"line",data:{labels:stages.map(s=>stageNames[s]),datasets:groups.map((g,i)=>({label:g,data:stages.map(s=>{const r=selected(g,s);return r?Number(r.median_patient_max_maf)*100:0;}),borderColor:i?"#526b4b":"#76543d",backgroundColor:i?"#879675":"#a67c59",tension:.25,pointRadius:4}))},options:{responsive:true,interaction:{mode:"index",intersect:false},plugins:{legend:{position:"bottom"},tooltip:{callbacks:{label:c:`${c.dataset.label}: ${fmt(c.raw)}%`}}},scales:{y:{beginAtZero:true,title:{display:true,text:"MAF máxima mediana (%)"}}}}});
 if(prevalenceChart)prevalenceChart.destroy();
 prevalenceChart=new Chart($("prevalence-chart"),{type:"bar",data:{labels:stages.map(s=>stageNames[s]),datasets:groups.map((g,i)=>({label:g,data:stages.map(s=>{const r=selected(g,s);return r?Number(r.prevalence)*100:0;}),backgroundColor:i?"#879675":"#a67c59",borderRadius:4}))},plugins:[{id:"wilsonIntervals",afterDatasetsDraw(chart){const ctx=chart.ctx;ctx.save();ctx.strokeStyle="#4f5549";ctx.lineWidth=1.4;chart.data.datasets.forEach((dataset,di)=>dataset.data.forEach((value,si)=>{const r=selected(groups[di],stages[si]);if(!r)return;const meta=chart.getDatasetMeta(di).data[si];if(!meta)return;const x=meta.x,lo=chart.scales.y.getPixelForValue(Number(r.prevalence_ci_low)*100),hi=chart.scales.y.getPixelForValue(Number(r.prevalence_ci_high)*100);ctx.beginPath();ctx.moveTo(x,hi);ctx.lineTo(x,lo);ctx.moveTo(x-4,hi);ctx.lineTo(x+4,hi);ctx.moveTo(x-4,lo);ctx.lineTo(x+4,lo);ctx.stroke();}));ctx.restore();}}],options:{responsive:true,interaction:{mode:"index",intersect:false},plugins:{legend:{position:"bottom"},tooltip:{callbacks:{label:c=>{const r=selected(groups[c.datasetIndex],stages[c.dataIndex]);return `${c.dataset.label}: ${fmt(c.raw,1)}% (IC95% ${fmt(Number(r?.prevalence_ci_low)*100,1)}–${fmt(Number(r?.prevalence_ci_high)*100,1)}%)`;}}}},scales:{y:{beginAtZero:true,max:100,title:{display:true,text:"Pacientes (%)"}}}}});
 const body=$("patient-tests-body");body.replaceChildren();
 (dashboard?.patient_tests||[]).forEach(r=>{const tr=document.createElement("tr");tr.innerHTML=`<td><b>${escapeHtml(r.Biological_group)}</b></td><td class="numeric">${fmt(r.n_patients,0)}</td><td class="numeric">${pFmt(r.friedman_p)}</td><td class="numeric">${pFmt(r.wilcoxon_p)}</td><td class="numeric">${fmt(r.presence_gained_AK_to_cSCC,0)} / ${fmt(r.presence_lost_AK_to_cSCC,0)}</td><td class="numeric">${pFmt(r.mcnemar_exact_p)}</td>`;body.append(tr);});
 if(!body.children.length){const tr=document.createElement("tr");tr.innerHTML='<td colspan="6">Métricas por paciente ainda não disponíveis.</td>';body.append(tr);}
}
function ensureTranscriptomicsSection(){
 if($("geo-transcriptomics"))return;
 const section=document.createElement("section");section.id="geo-transcriptomics";section.className="section geo-section";
 section.innerHTML=`<div class="section-kicker">04 / EVIDÊNCIA TRANSCRIPTÔMICA INDEPENDENTE</div>
 <div class="dashboard-heading"><div><h2>Da mutação à <em>expressão gênica.</em></h2><p>Coorte murina independente GSE158634 · RNA-seq bulk · 20 amostras depositadas.</p></div><span class="live-tag">GEO / MUS MUSCULUS</span></div>
 <div id="geo-status" class="status">Aguardando validação da matriz de contagens…</div>
 <div class="chart-grid"><article class="chart-card"><div class="chart-title"><div><span>CONTROLE DE QUALIDADE</span><h3>Profundidade das bibliotecas</h3></div></div><canvas id="geo-qc-chart" height="280"></canvas><p class="metric-explanation">Total de contagens por amostra (milhões). Marcadores sinalizados têm menos de 50% da mediana de tamanho de biblioteca do próprio estágio. São alertas para revisão, não exclusões automáticas.</p></article><article class="chart-card"><div class="chart-title"><div><span>EXPRESSÃO RELATIVA</span><h3>Genes do eixo ARID2–MAPK e Rora</h3></div></div><canvas id="geo-expression-chart" height="190"></canvas><p class="metric-explanation">Expressão resumida como log2(CPM + 0,5), com mediana por estágio. CPM corrige pelo tamanho da biblioteca; não é expressão absoluta nem substitui normalização/modelagem de RNA-seq para inferência.</p></article></div>
 <article class="chart-card"><div class="chart-title"><div><span>CONTRASTES DESCRITIVOS</span><h3>Diferença entre medianas de expressão</h3></div></div><div class="table-scroll"><table class="data-table"><thead><tr><th>Gene</th><th>Contraste</th><th>Δ mediana log2-CPM</th></tr></thead><tbody id="geo-contrasts-body"></tbody></table></div><p class="metric-explanation">Valores positivos indicam mediana maior no segundo estágio do contraste; negativos indicam menor. Não são log2 fold-changes modelados nem genes diferencialmente expressos.</p></article>
 <p class="data-note">A identificação individual dos animais não foi recuperada dos rótulos GEO. Por isso, esta camada não executa testes diferenciais, p-valores, enriquecimento funcional ou inferência sobre progressão individual. Os resultados humanos e murinos permanecem separados.</p>`;
 const anchor=document.querySelector(".team");(anchor||$("resultados")).before(section);
}
function drawTranscriptomics(){
 ensureTranscriptomicsSection();
 const status=dashboard?.geo_transcriptomics||{},rows=dashboard?.geo_expression||[];
 const label=status.status==="validated"?"Matriz validada":status.status==="counts_validated_gene_ids_unmapped"?"Contagens validadas; IDs gênicos sem mapeamento":status.status==="validation_failed"?"Validação não concluída":"Camada ainda não executada";
 const details=status.validation_message?": "+status.validation_message:"";
 $("geo-status").textContent=`${label} · ${status.n_samples??0}/20 amostras · ${status.n_features??0} features · ${status.n_target_genes??0} genes-alvo${details}`;
 $("geo-status").classList.toggle("error",status.status==="validation_failed");
 const qc=dashboard?.geo_sample_qc||[];\n if(geoQcChart)geoQcChart.destroy();\n geoQcChart=new Chart($("geo-qc-chart"),{type:"bar",data:{labels:qc.map(r=>r.sample_title),datasets:[{label:"Milhões de contagens",data:qc.map(r=>Number(r.library_size)/1e6),backgroundColor:qc.map(r=>r.qc_flag==="pass"?"#879675":"#b56b55"),borderRadius:3}]},options:{indexAxis:"y",responsive:true,plugins:{legend:{display:false},tooltip:{callbacks:{label:c=>`${fmt(c.raw,2)} milhões · ${qc[c.dataIndex]?.qc_flag==="pass"?"sem alerta":"revisar profundidade"}`}}},scales:{x:{beginAtZero:true,title:{display:true,text:"Milhões de contagens"}},y:{grid:{display:false}}}}});\n const genes=["Arid2","Rora","Egfr","Kras","Braf","Mapk1","Map2k1"];
 const available=genes.filter(g=>rows.some(r=>r.Gene===g));
 if(geoExpressionChart)geoExpressionChart.destroy();
 const stages=["Normal skin","Actinic keratosis","Cutaneous SCC"];
 geoExpressionChart=new Chart($("geo-expression-chart"),{type:"line",data:{labels:stages.map(s=>stageNames[s]||s),datasets:available.map((gene,i)=>({label:gene,data:stages.map(stage=>{const r=rows.find(x=>x.Gene===gene&&x.Stage===stage);return r?Number(r.median_log2_cpm):null;}),borderColor:["#76543d","#526b4b","#879675","#a67c59","#8c718d","#6f8991","#b18c52"][i%7],backgroundColor:"transparent",tension:.2,spanGaps:false}))},options:{responsive:true,interaction:{mode:"index",intersect:false},plugins:{legend:{position:"bottom"},tooltip:{callbacks:{label:c=>`${c.dataset.label}: ${fmt(c.raw,2)} log2-CPM`}}},scales:{y:{title:{display:true,text:"Mediana log2(CPM + 0,5)"}}}}});
 const body=$("geo-contrasts-body");body.replaceChildren();
 (dashboard?.geo_contrasts||[]).forEach(r=>{const tr=document.createElement("tr");tr.innerHTML=`<td><b>${escapeHtml(r.Gene)}</b></td><td>${escapeHtml(r.contrast)}</td><td class="numeric">${fmt(r.median_log2_cpm_difference,3)}</td>`;body.append(tr);});
 if(!body.children.length){const tr=document.createElement("tr");tr.innerHTML='<td colspan="3">Nenhum gene-alvo foi mapeado na matriz; não foi inferida correspondência por aproximação.</td>';body.append(tr);}
}
function renderExplorer(){drawMaf();}
function exportCsv(){
 const rows=filteredRows();if(!rows.length){showStatus("Não há registros para exportar com os filtros atuais.",true);return;}
 const headers=["Stage","Biological_group","Module","median"];
 const csv=[headers.join(";"),...rows.map(r=>headers.map(k=>String(r[k]??"").replace(/"/g,'""')).map(v=>`"${v}"`).join(";"))].join("\\r\\n");
 const blob=new Blob(["\\ufeff"+csv],{type:"text/csv;charset=utf-8"}),url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download="akstats-registros-filtrados.csv";a.click();URL.revokeObjectURL(url);
 showStatus(`Exportação preparada: ${rows.length} registro(s), respeitando os filtros ativos.`);
}
function setupInteractions(){
 const menu=$("menu-toggle"),nav=$("main-nav");if(menu&&nav){menu.addEventListener("click",()=>{const open=nav.classList.toggle("open");menu.setAttribute("aria-expanded",String(open));menu.textContent=open?"×":"☰";});nav.querySelectorAll("a").forEach(a=>a.addEventListener("click",()=>{nav.classList.remove("open");menu.setAttribute("aria-expanded","false");menu.textContent="☰";}));}
 ["stage-filter","module-filter"].forEach(id=>$(id).addEventListener("change",renderExplorer));$("gene-search").addEventListener("input",renderExplorer);
 $("reset-filters").addEventListener("click",()=>{$("stage-filter").value="all";$("module-filter").value="all";$("gene-search").value="";renderExplorer();});
 $("chart-toggle").addEventListener("click",()=>{chartMode=chartMode==="bar"?"line":"bar";$("chart-toggle").setAttribute("aria-pressed",String(chartMode==="line"));$("chart-toggle").textContent=chartMode==="line"?"Visualização em linha":"Visualização em barras";drawMaf();});
 $("export-csv").addEventListener("click",exportCsv);$("to-top").addEventListener("click",()=>window.scrollTo({top:0,behavior:"smooth"}));
 const observer=new IntersectionObserver(entries=>entries.forEach(e=>{if(e.isIntersecting){e.target.classList.add("is-visible");observer.unobserve(e.target);}}),{threshold:.12});
 document.querySelectorAll(".section,.stats-strip,.closing").forEach(el=>{el.classList.add("reveal");observer.observe(el);});
 window.addEventListener("scroll",()=> $("to-top").classList.toggle("visible",window.scrollY>600),{passive:true});
}
async function init(){
 setupInteractions();
 try{
  const response=await fetch(DATA_URL,{cache:"no-store"});if(!response.ok)throw new Error("Arquivo data/dashboard.json não encontrado. Execute o workflow de publicação dos dados.");
  dashboard=await response.json();fillKpis(dashboard.summary||{});
  [...new Set((dashboard.maf||[]).map(r=>r.Module).filter(Boolean))].forEach(m=>{const option=document.createElement("option");option.value=m;option.textContent=String(m).replace(/^\\d+\\. /,"");$("module-filter").append(option);});
  drawMaf();drawDelta();drawCooc();drawPatientMetrics();drawTranscriptomics();
  showStatus(`Atualizado em ${dashboard.generated_at?new Date(dashboard.generated_at).toLocaleString("pt-BR"):"data não informada"} · dados agregados do estudo publicado.`);
 }catch(err){showStatus(`Não foi possível carregar os dados: ${err.message}`,true);$("result-count").textContent="Dados indisponíveis";$("maf-table").innerHTML='<tr><td colspan="4">Os registros aparecerão assim que o workflow publicar o dashboard.json.</td></tr>';}
}
init();
