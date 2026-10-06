import React, { useEffect, useState } from 'react';
import { Link } from 'react-router-dom';
import './Home.css';

const tools = [
  {
    title: 'CholBindNet',
    eyebrow: 'Interpretable AI',
    description:
      'Predict and explore cholesterol-binding sites across protein structures with fast, interpretable results.',
    path: '/cholbindnet',
    img: '/CholThumb.webp',
    accent: 'cyan',
    featured: true,
  },
  {
    title: 'PDB Compare',
    eyebrow: 'Structure Analysis',
    description:
      'Compare protein structures and inspect residue-level distance changes in an interactive workflow.',
    path: '/pdbpaircompare',
    img: '/PDBCompare.webp',
    accent: 'violet',
  },
  {
    title: 'Current Flow Allostery',
    eyebrow: 'Allosteric Pathways',
    description:
      'Trace communication through protein structures and visualize the networks behind allosteric behavior.',
    path: '/newallosteric',
    img: '/Allosteric.webp',
    accent: 'lime',
  },
];

const ArrowIcon = () => (
  <svg aria-hidden="true" viewBox="0 0 24 24" fill="none">
    <path d="M5 12h14M13 6l6 6-6 6" stroke="currentColor" strokeWidth="1.8" />
  </svg>
);

let imagePreloadPromise;

function preloadToolImages() {
  if (!imagePreloadPromise) {
    imagePreloadPromise = Promise.all(
      tools.map(
        ({ img }) =>
          new Promise((resolve) => {
            const image = new Image();
            let settled = false;

            const finish = () => {
              if (settled) return;
              settled = true;

              if (typeof image.decode === 'function') {
                image.decode().catch(() => {}).finally(resolve);
              } else {
                resolve();
              }
            };

            image.onload = finish;
            image.onerror = finish;
            image.src = img;

            if (image.complete) finish();
          }),
      ),
    );
  }

  return imagePreloadPromise;
}

function Home() {
  const [imagesReady, setImagesReady] = useState(false);

  useEffect(() => {
    let isMounted = true;

    preloadToolImages().then(() => {
      if (isMounted) setImagesReady(true);
    });

    return () => {
      isMounted = false;
    };
  }, []);

  if (!imagesReady) {
    return (
      <main className="home-loading-screen" aria-label="Loading CompBio Tools">
        <span className="home-loading-dot" aria-hidden="true" />
      </main>
    );
  }

  return (
    <main className="home-shell">
      <div className="ambient ambient-one" />
      <div className="ambient ambient-two" />

      <nav className="site-nav" aria-label="Primary navigation">
        <Link className="brand" to="/" aria-label="CompBio Tools home">
          <span className="brand-mark" aria-hidden="true">
            <i />
            <i />
            <i />
          </span>
          <span>CompBio Tools</span>
        </Link>
      </nav>

      <section className="tools-hero" id="tools">
        <div className="intro-row">
          <div className="hero-copy">
            <p className="kicker"><span /> Computational biology, made interactive</p>
            <p className="hero-description">
              Research software for protein-structure analysis, allosteric pathways,
              and interpretable machine learning—built to turn complex data into
              useful scientific insight.
            </p>
          </div>

          <div className="intro-aside">
            <span className="availability">
              <span className="pulse" /> Three tools available
            </span>
            <p>
              Choose a workflow below. Each interface is purpose-built for
              computational biology research.
            </p>
          </div>
        </div>

        <div className="hero-graphic" aria-hidden="true">
          <div className="orbit orbit-one" />
          <div className="orbit orbit-two" />
          <div className="molecule-core">
            <span className="node node-a" />
            <span className="node node-b" />
            <span className="node node-c" />
            <span className="node node-d" />
            <span className="bond bond-a" />
            <span className="bond bond-b" />
            <span className="bond bond-c" />
          </div>
          <span className="data-chip chip-one">71.2% accuracy</span>
          <span className="data-chip chip-two">CPU inference</span>
        </div>

        <div className="tool-grid">
          {tools.map((tool) => (
            <Link
              to={tool.path}
              key={tool.path}
              className={`tool-card ${tool.featured ? 'featured' : ''}`}
            >
              <div className="card-image-wrap">
                <img
                  src={tool.img}
                  alt=""
                  loading="eager"
                  decoding="sync"
                  fetchPriority="high"
                />
                <div className="image-wash" />
                {tool.featured && <span className="featured-label">Featured</span>}
              </div>
              <div className="card-content">
                <p className={`card-eyebrow ${tool.accent}`}>{tool.eyebrow}</p>
                <h3>{tool.title}</h3>
                <p className="card-description">{tool.description}</p>
                <span className="card-cta">
                  Launch tool <ArrowIcon />
                </span>
              </div>
            </Link>
          ))}
        </div>
      </section>

      <footer>
        <p>Built for researchers working at the intersection of biology and computation.</p>
        <span>© {new Date().getFullYear()} CompBio Tools</span>
      </footer>
    </main>
  );
}

export default Home;
