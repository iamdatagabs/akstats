const menu=document.getElementById("menu-toggle"),nav=document.getElementById("main-nav");
if(menu&&nav){menu.addEventListener("click",()=>{const open=nav.classList.toggle("open");menu.setAttribute("aria-expanded",String(open));menu.textContent=open?"×":"☰";});nav.querySelectorAll("a").forEach(a=>a.addEventListener("click",()=>{nav.classList.remove("open");menu.setAttribute("aria-expanded","false");menu.textContent="☰";}));}
const topButton=document.getElementById("to-top");
if(topButton){topButton.addEventListener("click",()=>window.scrollTo({top:0,behavior:"smooth"}));window.addEventListener("scroll",()=>topButton.classList.toggle("visible",window.scrollY>550),{passive:true});}
const revealObserver=new IntersectionObserver(entries=>entries.forEach(entry=>{if(entry.isIntersecting){entry.target.classList.add("is-visible");revealObserver.unobserve(entry.target);}}),{threshold:.12});
document.querySelectorAll(".section,.question-section,.dashboard-teaser,.closing,.ticker").forEach(el=>{el.classList.add("reveal");revealObserver.observe(el);});
