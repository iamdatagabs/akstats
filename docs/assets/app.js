const DATA_URL = "data/dashboard.json";
let dashboard = null, mafChart = null, deltaChart = null;
const $ = (id) => document.getElementById(id);
const fmt = (n) => Number(n).toLocaleString("pt-BR", {maximumFractionDigits:2});
function showStatus(message, error=false){const el=$("data-status");el.textContent=message;el.style.color=error?"#9b4637":"";}
function fillKpis(s){$("kpi-genes").textContent=fmt(s.panel_genes??0);$("kpi-raw").textContent=fmt(s.source_variants??0);$("kpi-axis").textContent=fmt(s.axis_variants??0);$("kpi-patients").textContent=fmt(s.patients??0);}
function drawMaf(){
 const stage=$("stage-filter").value, group=$("module-filter").value;
 const rows=dashboard.maf.filter(r=>(stage==="all"||r.Stage===stage)&&(group==="all"||r.Module===group));
 const labels=[...new Set(rows.map(r=>r.Stage))];
 const datasets=["ARID2","MAPK"].map((g,i)=>({label:g,data:labels.map(s=>{const x=rows.find(r=>r.Stage===s&&r.Biological_group===g);return x?Number(x.median)*100:null;}),backgroundColor:i?"#879675":"#a67c59",borderRadius:5}));
 if(mafChart)mafChart.destroy();
 mafChart=new Chart($("maf-chart"),{type:"bar",data:{labels:labels.map(s=>({Normal:"Pele normal",AK:"Queratose actínica",SCC:"Carcinoma espinocelular"}[s]||s)),datasets},options:{responsive:true,plugins:{legend:{position:"bottom"}},scales:{y:{beginAtZero:true,title:{display:true,text:"MAF mediana (%)"}}}}});
}
function drawDelta(){
 const rows=dashboard.modules;
 if(deltaChart)deltaChart.destroy();
 deltaChart=new Chart($("delta-chart"),{type:"bar",data:{labels:rows.map(r=>r.Module.replace(/^\d+\. /,"")),datasets:[{label:"ΔMAF mediana (p.p.)",data:rows.map(r=>Number(r.median_delta_maf)*100),backgroundColor:"#879675",borderRadius:4}]},options:{indexAxis:"y",responsive:true,plugins:{legend:{display:false}},scales:{x:{title:{display:true,text:"Pontos percentuais"}}}}});
}
function drawCooc(){
 const host=$("cooc-summary");host.innerHTML="";
 dashboard.fisher.forEach(r=>{const row=document.createElement("div");row.className="cooc-row";const label=document.createElement("span");label.textContent=r.stage==="AK"?"Queratose actínica":"Carcinoma espinocelular";const value=document.createElement("b");value.textContent=`${r.both_present} pacientes com ARID2 + MAPK`;row.append(label,value);host.append(row);});
}
async function init(){
 try{
  const response=await fetch(DATA_URL,{cache:"no-store"});if(!response.ok)throw new Error("Arquivo de dados ainda não publicado.");
  dashboard=await response.json();fillKpis(dashboard.summary);
  const modules=[...new Set(dashboard.maf.map(r=>r.Module).filter(Boolean))];modules.forEach(m=>{const o=document.createElement("option");o.value=m;o.textContent=m.replace(/^\d+\. /,"");$("module-filter").append(o);});
  drawMaf();drawDelta();drawCooc();showStatus(`Atualizado em ${new Date(dashboard.generated_at).toLocaleString("pt-BR")} · fonte: suplementos do estudo publicado.`);
  $("stage-filter").addEventListener("change",drawMaf);
 }catch(err){showStatus(`Não foi possível carregar os dados: ${err.message} Execute o workflow de publicação.`,true);}
}
init();