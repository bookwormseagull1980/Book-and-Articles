"use strict";
const researchDiagram = document.getElementById("research-map");
if (researchDiagram) {
  const areas = [...document.querySelectorAll('map[name="research-navigation"] area')];
  const resizeMap = () => {
    if (!researchDiagram.naturalWidth) return;
    const scale = researchDiagram.clientWidth / researchDiagram.naturalWidth;
    for (const area of areas) {
      area.coords = area.dataset.coords.split(",").map(value => Math.round(Number(value) * scale)).join(",");
    }
  };
  researchDiagram.addEventListener("load", resizeMap);
  new ResizeObserver(resizeMap).observe(researchDiagram);
  resizeMap();
}
const search = document.getElementById("paper-search");
if (search) {
  const rows = [...document.querySelectorAll(".paper-row")];
  const status = document.getElementById("search-status");
  search.addEventListener("input", () => {
    const terms = search.value.toLowerCase().trim().split(/\s+/).filter(Boolean);
    let count = 0;
    for (const row of rows) {
      const haystack = row.dataset.search.toLowerCase();
      row.hidden = !terms.every(term => haystack.includes(term));
      if (!row.hidden) count += 1;
    }
    status.textContent = count === 0 ? "No matching papers" : count + (count === 1 ? " paper" : " papers");
  });
}
