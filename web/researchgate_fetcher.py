"""Fetch publication list from a ResearchGate profile page (HTML scraping)."""
from __future__ import annotations

import time

import requests
from bs4 import BeautifulSoup

USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)


def normalize_researchgate_url(url: str) -> str | None:
    """Extract or normalize ResearchGate profile URL."""
    url = (url or "").strip()
    if not url:
        return None
    if "researchgate.net" not in url.lower():
        return None
    if "/profile/" not in url.lower():
        return None
    if not url.startswith("http"):
        url = "https://" + url
    if "www." not in url and "researchgate.net" in url:
        url = url.replace("https://researchgate.net", "https://www.researchgate.net", 1)
    return url.split("?")[0].rstrip("/")


def get_researchgate_papers(
    profile_url: str,
    progress_callback: None = None,
) -> tuple[str, str, list[dict]]:
    """
    Fetch ResearchGate profile page and parse publications. Returns (profile_url, name, papers).
    papers = [{title, authors}]. May be incomplete if page is JS-rendered.
    """
    url = normalize_researchgate_url(profile_url)
    if not url:
        raise ValueError("Invalid ResearchGate profile URL")

    def report(msg: str) -> None:
        if progress_callback:
            progress_callback(msg)

    report("Fetching ResearchGate profile...")
    r = requests.get(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "text/html,application/xhtml+xml"},
        timeout=30,
    )
    r.raise_for_status()
    html = r.text
    soup = BeautifulSoup(html, "html.parser")

    name = "Unknown"
    title_tag = soup.find("title")
    if title_tag and title_tag.string:
        t = title_tag.string.strip()
        if " | " in t:
            name = t.split(" | ")[0].strip()
        else:
            name = t[:80]

    papers = []
    report("Parsing publications...")

    # Try multiple selector strategies (RG changes markup)
    pub_items = (
        soup.select("[data-testid='publication-item']")
        or soup.select(".nova-legacy-v-publication-item")
        or soup.select("div[class*='publication-item']")
        or soup.select("li[class*='publication']")
    )

    if not pub_items:
        # Fallback: any link to /publication/ with meaningful text
        for a in soup.select('a[href*="/publication/"]'):
            title = (a.get_text() or "").strip()
            if len(title) < 10 or len(title) > 500:
                continue
            papers.append({"title": title, "authors": []})
    else:
        for i, item in enumerate(pub_items):
            title_el = (
                item.select_one(".nova-legacy-v-publication-item__title")
                or item.select_one("[class*='publication-item__title']")
                or item.select_one("a[href*='/publication/']")
                or item.select_one("h2, h3, .title")
            )
            title = (title_el.get_text() if title_el else "").strip() if title_el else ""
            if not title and title_el and title_el.get("href"):
                title = (title_el.get_text() or "").strip()
            author_els = item.select(".nova-legacy-v-person-inline-item__fullname") or item.select(
                "[class*='person-inline'],[class*='author']"
            )
            authors = []
            for ae in author_els:
                a_text = (ae.get_text() or "").strip()
                if a_text and a_text not in authors:
                    authors.append(a_text)
            if title or authors:
                papers.append({"title": title or "(no title)", "authors": authors})
            if progress_callback and (i + 1) % 5 == 0:
                report(f"Parsed {i + 1}/{len(pub_items)} publications...")
            time.sleep(0.1)

    report(f"Done. Found {len(papers)} publications.")
    return url, name, papers
