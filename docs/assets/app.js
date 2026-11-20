const DATA_URL = "data/dashboard.json";
let dashboard = null, mafChart = null, deltaChart = null, chartMode = "bar";
const $ = (id) => document.getElementById(id);
const stageNames = {Normal:"Pele normal",AK:"Queratose actínica",SCC:"Carcinoma espinocelular"};
const fmt = (n, digits=2) => Number(n).toLocaleString("pt-BR",{maximumFractionDigits:digits});
const escapeHtml = (s) => String(s ?? "").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
function showStatus(message,error=false){const el=$("data-status");el.textContent=message;el.classList.toggle("error",error);}
function fillKpis(s){$("kpi-genes").textContent=fmt(s.panel_genes??0,0);$("kpi-raw").textContent=fmt(s.source_variants??0,0);$("kpi-axis").textContent=fmt(s.axis_variants??0,0);$("kpi-patients").textContent=fmt(s.patients??0,0);}
function filteredRows(){
 if(!dashboard?.maf)return [];
 const stage=$("stage-filter").value,module=$("module-filter").value,query=$("gene-search").value.trim().toLocaleLowerCase("pt-BR");
 return dashboard.maf.filter(r=>(stage==="all"||r.Stage===stage)&&(module==="all"||r.Module===module)&&(!query||[r.Module,r.Biological_group,r.Stage].join(" ").toLocaleLowerCase("pt-BR").includes(query)));
}
function drawMaf(){
 const rows=filteredRows(), stages=[...new Set(rows.map(r=>r.Stage))];
 const labels=stages.map(s=>stageNames[s]||s);
 const datasets=["ARID2","MAPK"].map((group,i)=>({label:group,data:stages.map(s=>{const item=rows.find(r=>r.Stage===s&&r.Biological_group===group);return item?Number(item.median)*100:null;}),backgroundColor:i?"#879675":"#a67c59",borderColor:i?"#526b4b":"#76543d",borderWidth:2,borderRadius:6,tension:.3,spanGaps:true}));
 if(mafChart)mafChart.destroy();
 mafChart=new Chart($("maf-chart"),{type:chartMode,data:{labels,datasets},options:{responsive:true,maintainAspectRatio:true,interaction:{mode:"index",intersect:false},plugins:{legend:{position:"bottom",labels:{usePointStyle:true,padding:20}},tooltip:{callbacks:{label:c=>`${c.dataset.label}: ${fmt(c.raw)}%`}}},scales:{y:{beginAtZero:true,title:{display:true,text:"MAF mediana (%)"},grid:{color:"#eee8de"}},x:{grid:{display:false}}},onClick:(event,elements)=>{if(!elements.length)return;const stage=stages[elements[0].index];$("stage-filter").value=stage;renderExplorer();}}});
 $("result-count").textContent=`${rows.length} registro(s) exibido(s) · ${new Set(rows.map(r=>r.Module)).size} módulo(s)`;
 renderTable(rows);
}
function drawDelta(){
 const rows=dashboard?.modules||[];
 if(deltaChart)deltaChart.destroy();
 deltaChart=new Chart($("delta-chart"),{type:"bar",data:{labels:rows.map(r=>String(r.Module||"").replace(/^\d+\. /,"")),datasets:[{label:"ΔMAF mediana (p.p.)",data:rows.map(r=>Number(r.median_delta_maf)*100),backgroundColor:"#879675",borderRadius:5}]},options:{indexAxis:"y",responsive:true,plugins:{legend:{display:false},tooltip:{callbacks:{label:c=>`${fmt(c.raw)} p.p.`}}},scales:{x:{title:{display:true,text:"Pontos percentuais"},grid:{color:"#eee8de"}},y:{grid:{display:false}}}}});
}
function drawCooc(){
 const host=$("cooc-summary");host.replaceChildren();
 (dashboard?.fisher||[]).forEach(r=>{const row=document.createElement("div");row.className="cooc-row";const label=document.createElement("span");label.textContent=stageNames[r.stage]||r.stage;const value=document.createElement("b");value.textContent=`${fmt(r.both_present,0)} pacientes com ARID2 + MAPK`;row.append(label,value);host.append(row);});
 if(!host.children.length)host.textContent="Sem registros de coocorrência disponíveis.";
}
function renderTable(rows){
 const body=$("maf-table");body.replaceChildren();
 if(!rows.length){const tr=document.createElement("tr");tr.innerHTML='<td colspan="4">Nenhum registro corresponde aos filtros atuais.</td>';body.append(tr);return;}
 rows.slice().sort((a,b)=>Number(b.median)-Number(a.median)).forEach(r=>{const tr=document.createElement("tr");tr.innerHTML=`<td>${escapeHtml(stageNames[r.Stage]||r.Stage)}</td><td><span class="group-pill ${r.Biological_group==="ARID2"?"arid":"mapk"}">${escapeHtml(r.Biological_group)}</span></td><td>${escapeHtml(String(r.Module||"—").replace(/^\d+\. /,""))}</td><td class="numeric">${fmt(Number(r.median)*100)}%</td>`;body.append(tr);});
}
function renderExplorer(){drawMaf();}
function exportCsv(){
 const rows=filteredRows();if(!rows.length){showStatus("Não há registros para exportar com os filtros atuais.",true);return;}
 const headers=["Stage","Biological_group","Module","median"];
 const csv=[headers.join(";"),...rows.map(r=>headers.map(k=>String(r[k]??"").replace(/"/g,'""')).map(v=>`"${v}"`).join(";"))].join("\r\n");
 const blob=new Blob(["\ufeff"+csv],{type:"text/csv;charset=utf-8"}),url=URL.createObjectURL(blob),a=document.createElement("a");
 a.href=url;a.download="akstats-registros-filtrados.csv";a.click();URL.revokeObjectURL(url);
 showStatus(`Exportação preparada: ${rows.length} registro(s), respeitando os filtros ativos.`);
}
function setupInteractions(){
 const menu=$("menu-toggle"),nav=$("main-nav");if(menu&&nav){menu.addEventListener("click",()=>{const open=nav.classList.toggle("open");menu.setAttribute("aria-expanded",String(open));menu.textContent=open?"×":"☰";});nav.querySelectorAll("a").forEach(a=>a.addEventListener("click",()=>{nav.classList.remove("open");menu.setAttribute("aria-expanded","false");menu.textContent="☰";}));}
 ["stage-filter","module-filter"].forEach(id=>$(id).addEventListener("change",renderExplorer));
 $("gene-search").addEventListener("input",renderExplorer);
 $("reset-filters").addEventListener("click",()=>{$("stage-filter").value="all";$("module-filter").value="all";$("gene-search").value="";renderExplorer();});
 $("chart-toggle").addEventListener("click",()=>{chartMode=chartMode==="bar"?"line":"bar";$("chart-toggle").setAttribute("aria-pressed",String(chartMode==="line"));$("chart-toggle").textContent=chartMode==="line"?"Visualização em linha":"Visualização em barras";drawMaf();});
 $("export-csv").addEventListener("click",exportCsv);
 $("to-top").addEventListener("click",()=>window.scrollTo({top:0,behavior:"smooth"}));
 const observer=new IntersectionObserver(entries=>entries.forEach(e=>{if(e.isIntersecting){e.target.classList.add("is-visible");observer.unobserve(e.target);}}),{threshold:.12});
 document.querySelectorAll(".section,.stats-strip,.closing").forEach(el=>{el.classList.add("reveal");observer.observe(el);});
 window.addEventListener("scroll",()=> $("to-top").classList.toggle("visible",window.scrollY>600),{passive:true});
}
async function init(){
 setupInteractions();
 try{
  const response=await fetch(DATA_URL,{cache:"no-store"});if(!response.ok)throw new Error("Arquivo data/dashboard.json não encontrado. Execute o workflow de publicação dos dados.");
  dashboard=await response.json();fillKpis(dashboard.summary||{});
  [...new Set((dashboard.maf||[]).map(r=>r.Module).filter(Boolean))].forEach(m=>{const option=document.createElement("option");option.value=m;option.textContent=String(m).replace(/^\d+\. /,"");$("module-filter").append(option);});
  drawMaf();drawDelta();drawCooc();
  showStatus(`Atualizado em ${dashboard.generated_at?new Date(dashboard.generated_at).toLocaleString("pt-BR"):"data não informada"} · fonte: suplementos do estudo publicado.`);
 }catch(err){showStatus(`Não foi possível carregar os dados: ${err.message}`,true);$("result-count").textContent="Dados indisponíveis";$("maf-table").innerHTML='<tr><td colspan="4">Os registros aparecerão assim que o workflow publicar o dashboard.json.</td></tr>';}
}
init();