"use strict";
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
