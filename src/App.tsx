import './App.css'

const stages = [
  'Oocyte collection',
  'Fertilization window',
  'Cleavage tracking',
  'Blastocyst scoring',
]

function App() {
  return (
    <main>
      <section className="hero">
        <nav className="topbar" aria-label="Primary">
          <a className="brand" href="/" aria-label="SURF 2026 Embryology home">
            <span className="brand-mark" aria-hidden="true" />
            <span>SURF 2026 Embryology</span>
          </a>
          <a className="nav-link" href="https://rishib.com">
            rishib.com
          </a>
        </nav>

        <div className="hero-grid">
          <div className="hero-copy">
            <p className="section-label">Caltech SURF 2026</p>
            <h1>Embryology research workspace for the 2026 SURF cycle.</h1>
            <p className="lede">
              A home for organizing experimental context, developmental stage
              references, imaging notes, and project updates as the embryology
              work takes shape.
            </p>
            <div className="actions">
              <a className="button primary" href="#timeline">
                View setup
              </a>
              <a className="button ghost" href="https://github.com/RB3572/surf2026Embryology">
                GitHub
              </a>
            </div>
          </div>

          <div className="visual-panel" aria-label="Embryology stage visual">
            <div className="microscope-field">
              <div className="cell cell-a" />
              <div className="cell cell-b" />
              <div className="cell cell-c" />
              <div className="cell cell-d" />
              <div className="scanline" />
            </div>
          </div>
        </div>
      </section>

      <section className="stage-band" id="timeline" aria-labelledby="timeline-title">
        <div>
          <p className="section-label">Project scaffold</p>
          <h2 id="timeline-title">Ready for research content</h2>
        </div>
        <div className="stage-list">
          {stages.map((stage, index) => (
            <article className="stage-card" key={stage}>
              <span>{String(index + 1).padStart(2, '0')}</span>
              <h3>{stage}</h3>
              <p>
                Placeholder structure for protocols, observations, figures, and
                notes once the research plan is finalized.
              </p>
            </article>
          ))}
        </div>
      </section>
    </main>
  )
}

export default App
